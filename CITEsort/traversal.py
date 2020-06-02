#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 20:58:29 2019

@author: lianqiuyu
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

#from Visualize import visualize_node,visualize_pair

class Traversal:
    
    def __init__(self,tree,c_type,method='bfs',nodelist=None,nodename=None,markers=None,n_samples=None,\
                 tree_summary=None,leaf_summary=None,n_components=None,ll=None,bic=None,leaf_ID=None,\
                 leaf_summary_code=None,multiplet_ratio=None,multiplet_predict=None):
        
        #print('initializing...')
        
        self.tree = tree
        self.method = method
        if self.method == 'bfs':
            self.nodelist = self.levelOrderTraversal()
        if self.method == 'dfs':
            self.nodelist = self.preorderTraversal()
        
        nodename_temp = ['_'.join(x.key) for x in self.nodelist]
        self.nodename = [str(i)+'_'+nodename_temp[i] for i in range(len(nodename_temp))]
        self.markers = [x[0] for x in self.nodelist[0].all_clustering[1]]    

        self.tree_summary, self.leaf_summary = self.summarize()
        
        self.n_components = self.leaf_summary.shape[0]
        self.ll = self.leaf_summary['ll'].sum()
        n_features = len(self.markers)
        
        mean_params = self.n_components * n_features
        if c_type == 'diag':
            cov_params = self.n_components * n_features
        if c_type == 'full':
            cov_params = self.n_components * n_features * (n_features + 1) / 2.
            
        n_parameters = int(self.n_components-1 + mean_params + cov_params)
        self.n_samples = len(self.nodelist[0].indices)
        self.bic = n_parameters * np.log(self.n_samples) - 2 * self.ll * self.n_samples
        self.leaf_ID = [int(x.split('_')[0]) for x in self.leaf_summary.index]
        self.predict_ACT_BCT()
        self.predict_multiplets()
        
        
        
    def summarize(self):
        #print('summarizing...')
        #num_node = len(self.nodename)
        n_samples = len(self.nodelist[0].indices)
        tree_summary = pd.DataFrame({'Count':[len(x.indices) for x in self.nodelist],
                                     'Proportion': [len(x.indices)/n_samples for x in self.nodelist],
                                     'Weight':[x.weight for x in self.nodelist],
                                     'll':[x.ll_tot for x in self.nodelist],
                                     'stop':[x.stop for x in self.nodelist]
                                     },index=self.nodename)
    
        mean_m = pd.DataFrame(np.zeros([tree_summary.shape[0],len(self.markers)]),
                              index = self.nodename,columns = self.markers)
        
        for i in range(mean_m.shape[0]):
            mean_m.iloc[i,:] = self.nodelist[i].mean_vec
        
        tree_summary = pd.concat([tree_summary,mean_m],axis=1)
        
        leaf_summary = tree_summary.loc[[x for x in self.nodename if x.split('_')[1]=='leaf'],:]
        leaf_summary = leaf_summary.sort_values(by='Count',ascending=False)
        return tree_summary,leaf_summary
    

    def get_node(self,nodeID):
        return self.nodelist[nodeID]
 
    
    def get_leaf_label(self):
        """generate label (one column, indicating which leaf cells are assigned.)"""
        label = pd.DataFrame({'GEM':self.tree.indices,'Label':[None]*len(self.tree.indices)},index=self.tree.indices)
        for i in range(len(self.nodename)):
            if self.nodename[i].split('_')[1] == 'leaf':
                label.loc[self.nodelist[i].indices,'Label'] = self.nodename[i]
                
        return label
    
    
    def plot_node(self,data,ID):
        node = self.nodelist[ID]
        node_data = data.loc[node.indices,:]
        plt.figure(figsize=(10,((data.shape[1]-1)//4+1)*2), dpi=96)
        plt.style.use('seaborn-white')
        if node.key == ('leaf',):
            for i in range(len(self.markers)):
                X = node_data.loc[:,self.markers[i]].values.reshape(-1, 1)
                bins = np.linspace(min(X),max(X),500)
                den = stats.norm.pdf(bins, node.mean_vec[i], np.sqrt(node.covariance_vec[i,i]))
                plt.subplot( (len(self.markers)-1)//5+1,5,i+1 )
                plt.hist(X,bins=30, density = True, color = "lightblue")
                plt.plot(bins,den,linewidth=1,color='black')
                plt.ylabel('density',fontsize=10)
                plt.title( self.markers[i],fontsize=12)
    
        else:
            
            for i in range(len(self.markers)):
                
                X = node_data.loc[:,self.markers[i]].values.reshape(-1, 1)
                bins = np.linspace(min(X),max(X),500)
                plt.subplot( (len(self.markers)-1)//4+1,5,i+1 )
                if (self.markers[i],) in node.all_clustering[1]:
                    weights = node.all_clustering[1][(self.markers[i],)]['component_weights']
                    means = node.all_clustering[1][(self.markers[i],)]['means']
                    covariances = node.all_clustering[1][(self.markers[i],)]['covariances']
                    y = np.zeros((len(bins),2))
                    y[:,0] = (weights[0] * stats.norm.pdf(bins, means[0], np.sqrt(covariances[0])))[:,0]
                    y[:,1] = (weights[1] * stats.norm.pdf(bins, means[1], np.sqrt(covariances[1])))[:,0]
                    if means[0] > means[1]:
                        cols = ['red','blue']
                    else:
                        cols = ['blue','red']
                    plt.plot(bins,y[:,0],linewidth=1,color=cols[0])
                    plt.plot(bins,y[:,1],linewidth=1,color=cols[1])
                else:
                    den = stats.norm.pdf(bins, node.mean_vec[i], np.sqrt(node.covariance_vec[i,i]))
                    plt.plot(bins,den,linewidth=1,color='black')
                
                plt.hist(X,bins=30, density = True, color = "lightblue")
                        
                subfig_title = self.markers[i] 
                if (self.markers[i],) == node.key:
                    plt.title( subfig_title,fontsize=12,color='red')
                else: 
                    plt.title( subfig_title,fontsize=12,color='darkgrey')
                    
                plt.ylabel('density',fontsize=10)
        
        plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9, hspace=0.4,wspace=0.45)
        plt.suptitle(self.nodename[ID]+' | '+str(len(node.indices))+' cells',fontsize=15,color="darkblue")
        plt.subplots_adjust(top=0.8)
        plt.show()   



    
    def predict_ACT_BCT(self):
        
        markers_cutoff = self._compute_markers_cutoff()
        markers = self.markers
        #leaf_p_markers = {}
        code = pd.DataFrame(np.zeros([len(self.leaf_ID),len(markers)]),index=self.leaf_summary.index,columns=markers)
        
        for leaf in code.index:
            mean = self.leaf_summary.loc[leaf,markers]
            for m in markers:
                code.loc[leaf,m] = 1 if mean[m] > markers_cutoff[m] else 0
        
        BCT_dic = {}
        ACT_dic = {}
        #ACT_tri_dic = {}
        
        for idx in code.index:
            #print(str(idx))
            if not BCT_dic:
                BCT_dic[idx] = np.sign(code.loc[idx,markers])
            else:
                new_center = np.sign(code.loc[idx,markers])
                new_flag = True
                terms = list(BCT_dic.keys())
                
                for i in range(len(terms)):
                    term1 = terms[i]
                    center1 = BCT_dic[term1]
                    for j in range(i,len(terms)):
                        term2 = terms[j]
                        center2 = BCT_dic[term2]
                        merge = pd.concat([center1,center2],axis=1)
                        if sum(new_center == merge.max(1)) == len(markers):
                            new_flag = False
                            if idx in ACT_dic:
                                ACT_dic[idx].append((term1,term2))
                            else:
                                ACT_dic[idx] = [(term1,term2)]
                
                if new_flag:
                    for i in range(len(terms)):
                        term1 = terms[i]
                        center1 = BCT_dic[term1]
                        for j in range(i+1,len(terms)):
                            term2 = terms[j]
                            center2 = BCT_dic[term2]
                            for k in range(j+1,len(terms)):
                                term3 = terms[k]
                                center3 = BCT_dic[term3]
                                merge = pd.concat([center1,center2,center3],axis=1)
                                if sum(new_center == merge.max(1)) == len(markers):
                                    new_flag = False
                                    if idx in ACT_dic:
                                        ACT_dic[idx].append((term1,term2,term3))
                                    else:
                                        ACT_dic[idx] = [(term1,term2,term3)]
                
                if new_flag:
                    BCT_dic[idx] = new_center
        
        
        leaf_summary_code = self.leaf_summary.drop(columns=self.markers)

        #leaf_summary_code.loc[list(ACT_dic.keys()),'Count'].sum()/data.shape[0]
        # 0.2504577309173559
        
        leaf_summary_code['BCT_predict'] = 0
        leaf_summary_code.loc[list(BCT_dic.keys()),'BCT_predict'] = 1
        
        leaf_summary_code['ACT_merge'] = None
        for term in ACT_dic.keys():
            leaf_summary_code.loc[term,'ACT_merge'] = str(ACT_dic[term])
        
        
        leaf_summary_code['merge_const'] = 0
        for term in ACT_dic.keys():
            p = []
            for pair in ACT_dic[term]:
                temp = 1
                for pair_i in pair:
                    temp = leaf_summary_code.loc[pair_i,'Weight'] * temp

                p.append(temp)
            #p /= len(ACT_dic[term]) 
            leaf_summary_code.loc[term,'merge_const'] = np.max(p)/leaf_summary_code.loc[term,'Weight']
        
        
        self.leaf_summary_code = pd.concat([leaf_summary_code,code],axis=1)


    
    def predict_multiplets(self):
        multiplet_predict = pd.Series([0]*self.n_samples,index=self.nodelist[0].indices)
        for leaf in self.leaf_summary_code.index:
            if self.leaf_summary_code.loc[leaf,'BCT_predict'] == 0 :
                multiplet_predict[self.nodelist[int(leaf.split('_')[0])].indices] = 1 
        
        self.multiplet_ratio = sum(multiplet_predict)/len(multiplet_predict)
        self.multiplet_predict = multiplet_predict
            
                
    def _compute_markers_cutoff(self):
        
        _all = self.nodelist[0]
        
        markers_cutoff = []
        for m in _all.all_clustering[1]:
            
            m1,m2 = _all.all_clustering[1][m]['means'][:,0]
            std1,std2 = np.sqrt(_all.all_clustering[1][m]['covariances'][:,0,0])
            s1,s2 = _all.all_clustering[1][m]['component_weights']
            inter_X = self._solve(m1,m2,std1,std2,s1,s2)
            if len(inter_X) == 1:
                markers_cutoff.append(inter_X)
            if (m1 - inter_X[0])*(m2 - inter_X[0]) < 0:
                markers_cutoff.append(inter_X[0])
            if (m1 - inter_X[1])*(m2 - inter_X[1]) < 0:
                markers_cutoff.append(inter_X[1])
        
        markers_cutoff = pd.Series(markers_cutoff,index=self.markers)
        
        return markers_cutoff
    
        
        
    def _solve(self,m1,m2,std1,std2,s1,s2):
        """solve equation: s1*N(m1,std1)=s2*N(m2,std2), return the intersection points of two weighted Gaussian"""
        a = 1/(2*std1**2) - 1/(2*std2**2)
        b = m2/(std2**2) - m1/(std1**2)
        c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log((std2*s1)/(std1*s2))
        return np.roots([a,b,c])



    # dfs
    def preorderTraversal(self):

        node = self.tree
        if node is None:
            return

        nodelist = []
        myStack = []

        while node or myStack:
            while node:
                nodelist.append(node)
                myStack.append(node)
                node = node.left
            node = myStack.pop()
            node = node.right   

        return nodelist


    # bfs
    def levelOrderTraversal(self): 
        #print('bfs...')
        node = self.tree
        if node is None: 
            return

        queue = [] 
        nodelist = []

        queue.append(node) 
        nodelist.append(node)

        while(len(queue) > 0): 
            node = queue.pop(0)         

            if node.left is not None: 
                nodelist.append(node.left)
                queue.append(node.left)

            if node.right is not None: 
                nodelist.append(node.right)
                queue.append(node.right) 

        return nodelist
















