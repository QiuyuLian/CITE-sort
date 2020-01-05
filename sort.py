#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 00:17:18 2019

@author: lianqiuyu
"""

from pomegranate import GeneralMixtureModel,MultivariateGaussianDistribution
import pandas as pd
import numpy as np
import itertools
from scipy import stats
from scipy.spatial import distance



def bhattacharyya(mu1, mu2, Sigma1, Sigma2):
    Sig = (Sigma1+Sigma2)/2
    ldet_s = np.linalg.det(Sig)
    ldet_s1 = np.linalg.det(Sigma1)
    ldet_s2 = np.linalg.det(Sigma2)
    d1 = distance.mahalanobis(mu1,mu2,np.linalg.inv(Sig))**2/8
    d2 = 0.5*np.log(ldet_s) - 0.25*np.log(ldet_s1) - 0.25*np.log(ldet_s2)
    return np.exp(-(d1+d2))



class BTree:

    def __init__(self, key, left = None, right = None, all_clustering = None,F_path=None,indices = None,\
                 prop=None,sample_weights=None,ll_tot=None,ll_vec=None,mean_vec=None,covariance_vec=None,\
                 ll_tot_delta=None,scores=None):
        self.key = key # a str
        self.right = right # a BstNode
        self.left = left # a BstNode
        self.indices = indices
        self.all_clustering = all_clustering
        self.F_path = F_path
        self.prop = prop
        self.sample_weights = sample_weights
        self.ll_tot = ll_tot
        self.ll_vec = ll_vec
        self.mean_vec = mean_vec
        self.covariance_vec = covariance_vec
        self.ll_tot_delta=ll_tot_delta
        self.scores = scores
            
    def display(self):
        lines, _, _, _ = self._display_aux()
        for line in lines:
            print(line)

    def _display_aux(self):
        """Returns list of strings, width, height, and horizontal coordinate of the root."""
        # No child.
        if self.right is None and self.left is None:
        #if self.right.key is 'leaf' and self.left.key is 'leaf':
            line = '%s' % '_'.join(self.key)
            width = len(line)
            height = 1
            middle = width // 2
            return [line], width, height, middle

        # Only left child.
        if self.right is None:
        #if self.right.key is 'leaf':
            lines, n, p, x = self.left._display_aux()
            s = '%s' % '_'.join(self.key)
            u = len(s)
            first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s
            second_line = x * ' ' + '/' + (n - x - 1 + u) * ' '
            shifted_lines = [line + u * ' ' for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, n + u // 2

        # Only right child.
        if self.left is None:
        #if self.left.key is 'leaf':
            lines, n, p, x = self.right._display_aux()
            s = '%s' % '_'.join(self.key)
            u = len(s)
            first_line = s + x * '_' + (n - x) * ' '
            second_line = (u + x) * ' ' + '\\' + (n - x - 1) * ' '
            shifted_lines = [u * ' ' + line for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, u // 2

        # Two children.
        left, n, p, x = self.left._display_aux()
        right, m, q, y = self.right._display_aux()
        s = '%s' % '_'.join(self.key)
        u = len(s)
        first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s + y * '_' + (m - y) * ' '
        second_line = x * ' ' + '/' + (n - x - 1 + u + y) * ' ' + '\\' + (m - y - 1) * ' '
        if p < q:
            left += [n * ' '] * (q - p)
        elif q < p:
            right += [m * ' '] * (p - q)
        zipped_lines = zip(left, right)
        lines = [first_line, second_line] + [a + u * ' ' + b for a, b in zipped_lines]
        return lines, n + m + u, max(p, q) + 2, n + u // 2



def ScanFeatures(data,F_set,root,fix_k,ndim):
    
    ########################################
    # function:
    #       generate separable feature candidates in ndim-D space
    #
    ########################################
    
    separability = {1:0.1,2:0.5}
    all_clustering = {}
    separable_features = []
    
    for item in itertools.combinations(F_set, ndim):
        
        if item in root.F_path:
            continue
        
        x = data.loc[:,item]
        subgmm = GeneralMixtureModel.from_samples(MultivariateGaussianDistribution, n_components=fix_k, X=x,weights=w)
        
        sub_clustering = {}
        sub_clustering['means'] = np.array([subgmm.distributions[k].parameters[0] for k in range(fix_k)])
        sub_clustering['covariances'] = np.array([subgmm.distributions[k].parameters[1] for k in range(fix_k)])
        sub_clustering['component_weights'] = np.exp(subgmm.weights)
        
        sub_clustering['bhat_similarity'] = bhattacharyya(sub_clustering['means'][0,:],sub_clustering['means'][1,:],sub_clustering['covariances'][0,:],sub_clustering['covariances'][1,:])
        
        #sub_clustering['entropy'] = stats.entropy(sub_clustering['component_weights'])
        #sub_clustering['posterior'] = subgmm.predict_proba(x)
        
        all_clustering[item] = sub_clustering
            
        if sub_clustering['bhat_similarity'] < separability[ndim]:
            if item not in root.F_path:
                separable_features.append(item)
                sub_clustering['posterior'] = subgmm.predict_proba(x)
                    
    root.all_clustering[ndim] = all_clustering
    
    return separable_features
    


def ScoreFeatures(data,root,fix_k,max_ndim,rescan_cut):
    
    separable_features = []
    
    w = root.sample_weights
    root.mean_vec = (data * w.values.reshape(-1,1)).sum() / w.sum()
    root.covariance_vec = (data**2 * w.values.reshape(-1,1)).sum() / w.sum() - root.mean_vec**2
    
    ll_vec = loglikelihood(data,w,prop, mean_vec = root.mean_vec, covariance_vec = root.covariance_vec)   ####
    root.ll_tot = np.sum(ll_vec) 
    root.ll_vec = pd.Series(ll_vec, index=data.columns).sort_values(ascending=True)
    root.all_clustering = {}
    
    ndim = 1
    F_set = data.columns.values.tolist()
    separable_features = ScanFeatures(data,F_set,root,fix_k,ndim)
    
    if len(separable_features)==0:
        F_set = []
        for item in root.all_clustering[ndim]:
            val = root.all_clustering[ndim][item]['bhat_similarity']
            if val < 0.5:
                F_set.append(item[0])

        if len(F_set) > 0:
            for ndim in range(2,max_ndim+1):
                separable_features = ScanFeatures(data,F_set,root,fix_k,ndim)
                if len(separable_features) > 0:
                    break
                
        
    if len(separable_features) > 0:
        
        scores = pd.Series(np.zeros(len(separable_features)),index=separable_features)
        
        for f in separable_features:
            #print(f)
            posterior = root.all_clustering[ndim][f]['posterior']
            component_weights = root.all_clustering[ndim][f]['component_weights']
            ll_vec1 = loglikelihood(data, posterior[:,0]* w, prop*component_weights[0])
            ll_vec2 = loglikelihood(data, posterior[:,1]* w, prop*component_weights[1])
            scores[f] = sum(ll_vec1) + sum(ll_vec2) - root.ll_tot 
        
        root.scores = scores
    
    return separable_features



def loglikelihood(data,w,prop,mean_vec=None,covariance_vec=None):
    
    if mean_vec is None:
        mean_vec = (data * w.values.reshape(-1,1)).sum() / w.sum()
        covariance_vec = (data**2 * w.values.reshape(-1,1)).sum() / w.sum() - mean_vec**2
    
    ll_vec = prop * np.array(list(map(lambda i: np.sum(w*stats.norm.logpdf(data.iloc[:,i],mean_vec[i],np.sqrt(covariance_vec[i])))/data.shape[0],range(data.shape[1]))))
    return ll_vec
    


def sort(data,w,F_path=[],prop=1,rescan_cut=0.7,tol=1e-2,fix_k=2,max_ndim=2):
    
    ## construct a node to save basic info about current conditional probability.
    root = BTree('tmp')
    root.prop = prop
    root.F_path = F_path # 
    root.sample_weights = w
        
    separable_features = ScoreFeatures(data,root,fix_k,max_ndim,rescan_cut)

    
    if len(separable_features)==0:
        root.key = ('leaf',)
        return root
    
    best_feature = root.scores.idxmax()
    root.key = best_feature
    root.ll_tot_delta = root.scores.max()
    
    if root.scores.max() < tol:
        root.key = ('leaf',)
        return root
    
    ndim = len(best_feature)
    posterior = root.all_clustering[ndim][best_feature]['posterior']
    component_weights = root.all_clustering[ndim][best_feature]['component_weights']
    
    ## branch cells, component with higher mean goes right.
    p0_mean = root.all_clustering[ndim][best_feature]['means'][0,:]
    p1_mean = root.all_clustering[ndim][best_feature]['means'][1,:]
    
    flag = True
    
    if ndim == 1:
        flag = p0_mean[0] < p1_mean[0]
    else:
        p0_cosine = sum(p0_mean)/np.sqrt(sum(p0_mean**2))
        p1_cosine = sum(p1_mean)/np.sqrt(sum(p1_mean**2))
        flag = p0_cosine < p1_cosine

    if flag:
        root.left = sort(data,posterior[:,0] * w,F_path + [best_feature],component_weights[0]*prop,rescan_cut, tol,fix_k,max_ndim)
        root.right = sort(data,posterior[:,1] * w,F_path + [best_feature],component_weights[1]*prop, rescan_cut,tol,fix_k,max_ndim)

    else:
        root.left = sort(data,posterior[:,1] * w,F_path + [best_feature],component_weights[1]*prop,rescan_cut, tol,fix_k,max_ndim)
        root.right = sort(data,posterior[:,0] * w,F_path + [best_feature],component_weights[0]*prop,rescan_cut, tol,fix_k,max_ndim)

    return root




data = pd.read_csv('./data/cellhashing_KC.csv',header=0,index_col=0)

#data_2 = data**2

#F_set = data.columns.values.tolist()

import random


random.seed( 30000 )




rescan_cut=0.8
fix_k = 2
ndim = 1
max_ndim = 1
prop = 1
tol = 1e-2
F_path = []
w = pd.Series(np.array([1]*data.shape[0]),index=data.index)


import time



start_time = time.time()
tree = sort(data[['CD3','CD19','CD8a','CD11c','CD14','CD16','CD56','CD4']],w,F_path,prop,rescan_cut,tol=1e-2)
print("--- %s seconds ---" % (time.time() - start_time))
tree.display()





# sort(data,w,F_path=[],prop=1,rescan_cut=0.7,tol=1e-2,fix_k=2,max_ndim=2)



from traversal import Traversal
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.mixture import GaussianMixture

multiplets = np.genfromtxt('./data/multiplets.csv',dtype='str')
msm_labels = pd.Series([0]*data.shape[0],index=data.index)
msm_labels[multiplets] = 1




tree_tr = Traversal(tree)
print(tree_tr.ll_tot)


labels = tree_tr.leaf_label
posteriors = tree_tr.posterior
props = tree_tr.prop

len(np.unique(labels))
n_cluster = len(np.unique(labels))

print(adjusted_mutual_info_score(labels, msm_labels))



label_name = labels.unique()
labels_msm_df = pd.DataFrame(np.zeros([len(label_name),2]),columns=['counts','msm'],index=label_name)

for i in range(len(label_name)):
    idx = labels == label_name[i]
    labels_msm_df.iloc[i,0] = sum(idx)
    labels_msm_df.iloc[i,1] = sum(msm_labels[idx]==1) 
    

labels_msm_df['msm_ratio'] = labels_msm_df['msm']/labels_msm_df['counts']
labels_msm_df['ratio'] = labels_msm_df['counts']/data.shape[0]


labels_msm_df_sort = labels_msm_df.sort_values(by=['msm_ratio'],ascending=False)


centers_df = pd.DataFrame(np.zeros([len(label_name),10]),index=label_name,columns=data.columns)

for leaf in label_name:
    idx = tree_tr.nodename.index(leaf)
    centers_df.loc[leaf,:] = tree_tr.nodelist[idx].mean_vec


stat = pd.concat([labels_msm_df,centers_df], axis=1, sort=False)
stat_sort = stat.sort_values(by=['counts'],ascending=False)




data_tsne = pd.read_csv('./result/kc_complete/tsne/data_tsne.csv',header=0,index_col=0)

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib import cm
import seaborn as sns

def showLabel(data,labels_pred):

    labels_set = np.unique(labels_pred)
    k = len(labels_set)

    #cmap = cm.get_cmap('Set3')

    cell_colors_dic = dict(zip(labels_set,[cm.rainbow(i) for i in np.linspace(0, 1, k)]))
    #cell_colors = pd.DataFrame({'label':[None]*data.shape[0]},index=data.index)
    sns.set(style="white")

    fig, ax = plt.subplots()
    for i in labels_set:
        plt.plot(data[labels_pred==i,0],data[labels_pred==i,1],'o',
                 c=cell_colors_dic[i],label=i,markersize=0.6,alpha=1)
    #plt.xlabel('M0',fontsize=12)
    #plt.ylabel('M1',fontsize=12)
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.axis('tight')
    plt.title('number of clusters: %d' % k,fontsize=15)

    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,markerscale=8)   
    plt.show()
    
    

showLabel(data_tsne.values,labels)


target = labels=='67_leaf'
target_labels = pd.Series([0]*data.shape[0],index=data.index)
target_labels[target] = 1

showLabel(data_tsne.values,target_labels)





start_time = time.time()
gmm = GaussianMixture(n_cluster).fit(data)
print("--- %s seconds ---" % (time.time() - start_time))

gmm_lables = gmm.predict(data)


print(gmm.lower_bound_)
print(adjusted_mutual_info_score(gmm_lables, msm_labels))







gmm_lables_s = pd.Series(gmm_lables,index=data.index)
gmm_label_name = gmm_lables_s.unique()

gmm_labels_msm_df = pd.DataFrame(np.zeros([len(gmm_label_name),2]),columns=['counts','msm'],index=gmm_label_name)

for i in range(len(gmm_label_name)):
    idx = gmm_lables_s == gmm_label_name[i]
    gmm_labels_msm_df.iloc[i,0] = sum(idx)
    gmm_labels_msm_df.iloc[i,1] = sum(msm_labels[idx]==1) 

gmm_labels_msm_df['msm_ratio'] = gmm_labels_msm_df['msm']/gmm_labels_msm_df['counts']
gmm_labels_msm_df['ratio'] = gmm_labels_msm_df['counts']/data.shape[0]
gmm_labels_msm_df_sort = gmm_labels_msm_df.sort_values(by=['msm_ratio'],ascending=False)






record_df_r = pd.DataFrame({'method':['RGMM']*n_cluster+['GMM']*n_cluster,
                            'value':labels_msm_df_sort['msm_ratio'].values.tolist()+gmm_labels_msm_df_sort['msm_ratio'].values.tolist()})

plt.figure(figsize=(3,3), dpi=96)
ax = sns.violinplot(x="method",y="value",data=record_df_r,showfliers=False,color='grey',width=0.5,linewidth=1)
ax = sns.swarmplot(x="method",y="value",data=record_df_r,color='black',dodge=True,size=4)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylabel('msm ratio')
plt.yticks([0,1])
#plt.savefig('ADT_RNA_cor.pdf') 
plt.show()




