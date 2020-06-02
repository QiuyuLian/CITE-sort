
import sys
sys.path.append("./CITEsort")
import numpy as np
from sklearn.mixture import GaussianMixture
import itertools
from scipy import stats
import operator
from scipy.spatial import distance
from BTree import BTree
import copy
#from scipy.signal import upfirdn
#import pandas as pd


def Matryoshka(data,merge_cutoff=0.1,max_k=10,max_ndim=2,bic='bic'):
    
    if data.shape[0] <= 20:
        root = BTree(('leaf',))
        root.indices = data.index.values.tolist()
        root.all_clustering_dic = _set_small_leaf(data)
        return root
    
    separable_features, bipartitions, scores, all_clustering_dic = HiScanFeatures(data,merge_cutoff,max_k,max_ndim,bic)
    
    if len(separable_features) == 0:
        root = BTree(('leaf',))
        root.indices = data.index.values.tolist()
        root.all_clustering_dic = all_clustering_dic
        return root
    
    idx_best = np.argmax(scores)
    best_feature = separable_features[idx_best]
    best_partition = bipartitions[best_feature]
    
    ## construct current node  
    root = BTree(best_feature)
    root.indices = data.index.values.tolist()
    root.all_clustering_dic = all_clustering_dic
    #root.marker_summary = marker_summary
    #root.para = para

    ## branch cells, component with higher mean goes right.
    p1_mean = data.loc[best_partition, best_feature].mean(0)
    p2_mean = data.loc[~best_partition, best_feature].mean(0)
    
    flag = True
    if len(p1_mean) == 1:
        flag = p1_mean.values > p2_mean.values
    else:
        p1_cosine = sum(p1_mean)/np.sqrt(sum(p1_mean**2))
        p2_cosine = sum(p2_mean)/np.sqrt(sum(p2_mean**2))
        flag = p1_cosine > p2_cosine

    if flag:
        child_right = data.iloc[best_partition, :]
        child_left = data.iloc[~best_partition, :] 
        root.where_dominant = 'right'
    else:
        child_right = data.iloc[~best_partition, :]
        child_left = data.iloc[best_partition, :]
        root.where_dominant = 'left'
    
    ## recursion
    root.left = Matryoshka(child_left,merge_cutoff,max_k,max_ndim,bic)
    root.right = Matryoshka(child_right,merge_cutoff,max_k,max_ndim,bic)

    return root



def HiScanFeatures(data,merge_cutoff,max_k,max_ndim,bic):
    
    ndim = 1
    all_clustering_dic = {}
    separable_features, bipartitions, scores, all_clustering_dic[ndim] = ScanFeatures(data,merge_cutoff,max_k,ndim,bic)
    
    if len(separable_features) == 0:

        rescan_features = []
        for item in all_clustering_dic[ndim]:
            val = all_clustering_dic[ndim][item]['similarity_stopped']
            if val > 0.1 and val < 0.5:
                rescan_features.append(item[0])
                
        for ndim in range(2,max_ndim+1):
            separable_features, bipartitions, scores, all_clustering_dic[ndim] = ScanFeatures(data[rescan_features],merge_cutoff,max_k,ndim,bic)
            if len(separable_features) >= 1:
                break
        
    return separable_features, bipartitions, scores, all_clustering_dic
    


def ScanFeatures(data,merge_cutoff,max_k,ndim,bic):
    
    F_set = data.columns.values.tolist()
    
    all_clustering = {}
    separable_features = []
    bipartitions = {}
    scores = []
        
    for item in itertools.combinations(F_set, ndim):
        x = data.loc[:,item]
        all_clustering[item] = Clustering(x,merge_cutoff,max_k,bic)
    
    for item in all_clustering:
        if all_clustering[item]['mp_ncluster'] > 1:
            separable_features.append(item)
            bipartitions[item] = all_clustering[item]['max_ent_p']
            scores.append(all_clustering[item]['max_ent'])
            
    return separable_features, bipartitions, scores, all_clustering



def Clustering(x,merge_cutoff,max_k,bic):
    
    val,cnt = np.unique(x.values.tolist(),return_counts=True)
    
    if len(val) < 50:
        clustering = _set_one_component(x) 
   
    else:
    
        k_bic,_ = BIC(x,max_k,bic)
    
        if k_bic == 1:    
            # if only one component, set values
            clustering = _set_one_component(x)      
        else:
            
            bp_gmm = GaussianMixture(k_bic).fit(x)
            clustering = merge_bhat(x,bp_gmm,merge_cutoff)
            
            if clustering['mp_ncluster'] > 1:
    
                merged_label = clustering['mp_clustering']
                labels, counts = np.unique(merged_label, return_counts=True)
                per = counts/np.sum(counts)
                ents = [stats.entropy([per_i, 1-per_i],base=2) for per_i in per]
                clustering['max_ent'] = np.max(ents)
                best_cc_idx = np.argmax(ents)
                best_cc_label = labels[best_cc_idx]
                clustering['max_ent_p'] = merged_label == best_cc_label
        
    return clustering



def bhattacharyya_dist(mu1, mu2, Sigma1, Sigma2):
    Sig = (Sigma1+Sigma2)/2
    ldet_s = np.linalg.det(Sig)
    ldet_s1 = np.linalg.det(Sigma1)
    ldet_s2 = np.linalg.det(Sigma2)
    d1 = distance.mahalanobis(mu1,mu2,np.linalg.inv(Sig))**2/8
    d2 = 0.5*np.log(ldet_s) - 0.25*np.log(ldet_s1) - 0.25*np.log(ldet_s2)
    return d1+d2



def merge_bhat(x,bp_gmm,cutoff):

    clustering = {}
    clustering['bp_ncluster'] = bp_gmm.n_components
    clustering['bp_clustering'] = bp_gmm.predict(x)
    clustering['bp_pro'] = bp_gmm.weights_
    clustering['bp_mean'] = bp_gmm.means_
    clustering['bp_Sigma'] = bp_gmm.covariances_
    
    #clustering['last_pair_similarity'] = _get_last_pair_similarity_2D(x,bp_gmm)
    gmm = copy.deepcopy(bp_gmm) 
    
    mu = gmm.means_
    Sigma = gmm.covariances_
    weights = list(gmm.weights_)
    posterior = gmm.predict_proba(x)
    
    current_ncluster = len(mu)
    mergedtonumbers = [int(item) for item in range(current_ncluster)]

    merge_flag = True
    clustering['bhat_dic_track'] = {}
    merge_time = 0

    while current_ncluster > 1 and merge_flag:

        bhat_dic = {}

        for c_pair in itertools.combinations(range(current_ncluster), 2):
            m1 = mu[c_pair[0],:]
            m2 = mu[c_pair[1],:]
            Sigma1 = Sigma[c_pair[0],:,:]
            Sigma2 = Sigma[c_pair[1],:,:]
            bhat_dic[c_pair] = np.exp(-bhattacharyya_dist(m1, m2, Sigma1, Sigma2))

        clustering['bhat_dic_track'][merge_time] = bhat_dic
        merge_time = merge_time + 1
        
        max_pair = max(bhat_dic.items(), key=operator.itemgetter(1))[0]
        max_val = bhat_dic[max_pair]

        if max_val > cutoff:
            merged_i,merged_j = max_pair
            # update mergedtonumbers
            for idx,val in enumerate(mergedtonumbers):
                if val == merged_j:
                    mergedtonumbers[idx] = merged_i
                if val > merged_j:
                    mergedtonumbers[idx] = val - 1
                    
            # update parameters
            weights[merged_i] = weights[merged_i] + weights[merged_j]
            
            posterior[:,merged_i] = posterior[:,merged_i] + posterior[:,merged_j]
            
            w = posterior[:,merged_i]/np.sum(posterior[:,merged_i])
            mu[merged_i,:] = np.dot(w,x)# update                                 
            
            x_centered = x.apply(lambda xx: xx-mu[merged_i,:],1)
            Sigma[merged_i,:,:] = np.cov(x_centered.T,aweights=w,bias=1)

            del weights[merged_j]
            #weights = np.delete(weights,merged_j,0)
            mu = np.delete(mu,merged_j,0)
            Sigma = np.delete(Sigma,merged_j,0)
            posterior = np.delete(posterior,merged_j,1)
            current_ncluster = current_ncluster - 1

        else:
            merge_flag = False
    
    
    clustering['similarity_stopped'] = np.min(list(bhat_dic.values()))
    clustering['mp_ncluster'] = mu.shape[0]
    clustering['mergedtonumbers'] = mergedtonumbers
    clustering['mp_clustering'] = list(np.apply_along_axis(np.argmax,1,posterior))
    
    return clustering



def _set_small_leaf(data):
    all_clustering_dic = {}
    all_clustering_dic[1] = {}
    
    F_set = data.columns.values.tolist()
    all_clustering = {}

    for item in itertools.combinations(F_set, 1):
        x = data.loc[:,item]
        all_clustering[item] = _set_one_component(x)
    
    all_clustering_dic[1] = all_clustering
    
    return all_clustering_dic



def _set_one_component(x):
    
    clustering = {}
    clustering['bp_ncluster'] = 1
    clustering['bp_clustering'] = [0]*len(x)
    clustering['bp_pro'] = [1]
    clustering['bp_mean'] = np.mean(x)
    clustering['bp_Sigma'] = np.var(x)
    clustering['bhat_dic_track'] = {}
    clustering['similarity_stopped'] = 1
    clustering['mp_ncluster'] = 1
    clustering['mp_clustering'] = [0]*len(x)
    clustering['mergedtonumbers'] = [0]
    
    return clustering



def BIC(X, max_k = 10,bic = 'bic'):
    """return best k chosen with BIC method"""
    
    bic_list = _get_BIC_k(X, min(max_k,len(np.unique(X))))
    
    if bic == 'bic':   
        return min(np.argmin(bic_list)+1,_FindElbow(bic_list)),bic_list
    elif bic == 'bic_min':   
        return np.argmin(bic_list)+1,bic_list
    elif bic == 'bic_elbow':
        return _FindElbow(bic_list),bic_list

    
    
def _get_BIC_k(X, max_k):
    """compute BIC scores with k belongs to [1,max_k]"""
    bic_list = []
    for i in range(1,max_k+1):
        gmm_i = GaussianMixture(i).fit(X)
        bic_list.append(gmm_i.bic(X))
    return bic_list



def _FindElbow(bic_list):
    """return elbow point, defined as the farthest point from the line through the first and last points"""
    if len(bic_list) == 1:
        return 1
    else:
        a = bic_list[0] - bic_list[-1]
        b = len(bic_list) - 1
        c = bic_list[-1]*1 - bic_list[0]*len(bic_list)
        dis = np.abs(a*range(1,len(bic_list)+1) + b*np.array(bic_list) + c)/np.sqrt(a**2+b**2)
        return np.argmax(dis)+1


