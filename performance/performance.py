#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 22:39:00 2020

@author: lianqiuyu
"""

import seaborn as sns
from matplotlib import pyplot as plt
#from sort_harddivision import sort
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture,BayesianGaussianMixture


sns.set(style="whitegrid")
import time
from ReSplit import ReSplit
from BTreeTraversal import BTreeTraversal
#from DEmerge import DEmerge

namelist = ['PBMC_1k','PBMC_1k_b','PBMC_2k', 'PBMC_5k', 'PBMC_8k', 'MALT_8k', 'CBMC_8k','PBMC_16k']
datapath = './datasets/data_for_performance'
savepath = './performance'



from sys import argv



max_cluster_num = 50

def find_k(data,c_type,max_cluster_num=100):
    k_list = []
    inertia = []
    
    for k in range(1, max_cluster_num + 1):
    
        gmm = GaussianMixture(k,covariance_type=c_type).fit(data)
        k_list.append(k)
        inertia.append(gmm.bic(data))
    
    idx = np.argmin(inertia)
    final_k = k_list[idx]
    return final_k


merge_cutoff = 0.1
record_full ={}

#for i in range(len(namelist)):
    
#name = namelist[i]
name = argv[1]
print(name)
data = pd.read_csv(datapath+'/'+name+'_ADT_clr_10markers.csv',header=0,index_col=0)
#N=data.shape[0]

record_sort = pd.DataFrame(np.zeros([10,3]),columns = ['time','ll','n_component'])
record_gmm = pd.DataFrame(np.zeros([10,3]),columns = ['time','ll','n_component'])
record_ngmm = pd.DataFrame(np.zeros([10,3]),columns = ['time','ll','n_component'])
record_dpgmm = pd.DataFrame(np.zeros([10,3]),columns = ['time','ll','n_component'])

record_gmm_fix_k = pd.DataFrame(np.zeros([10,3]),columns = ['time','ll','n_component'])
record_ngmm_fix_k = pd.DataFrame(np.zeros([10,3]),columns = ['time','ll','n_component'])



for t in range(10):
    
    print('CITE-sort.')
    start_time = time.time()
    rgmm = ReSplit(data,merge_cutoff) #sort(data,F_path,N,c_type,weight,rescan_cut,tol,fix_k,max_ndim)
    record_sort.iloc[t,0]  = time.time() - start_time
    # print("--- %s seconds ---" % (t1))
    trav = BTreeTraversal(rgmm)
    record_sort.iloc[t,1]  = trav.ll
    #record_sort.iloc[t,2]  = trav.bic
    record_sort.iloc[t,2] =  trav.n_components
    #print(rgmm_ll)
    
    print('full GMM.')
    start_time = time.time()
    final_k = find_k(data,'full', max_cluster_num)
    gmm = GaussianMixture(final_k).fit(data)
    record_gmm.iloc[t,0] = time.time() - start_time
    #print("--- %s seconds ---" % gmm_time)
    #record_gmm.iloc[t,2] = gmm.bic(data)
    record_gmm.iloc[t,1] = gmm.score(data)
    record_gmm.iloc[t,2] = final_k #trav.n_components
    
    
    print('full GMM wigh fix k.')
    start_time = time.time()
    gmm = GaussianMixture(trav.n_components).fit(data)
    record_gmm_fix_k.iloc[t,0] = time.time() - start_time
    #print("--- %s seconds ---" % gmm_time)
    #record_gmm.iloc[t,2] = gmm.bic(data)
    record_gmm_fix_k.iloc[t,1] = gmm.score(data)
    record_gmm_fix_k.iloc[t,2] = trav.n_components
    

    print('naive GMM.')
    start_time = time.time()
    final_k = find_k(data,'diag',max_cluster_num)
    ngmm = GaussianMixture(final_k,covariance_type='diag').fit(data)
    record_ngmm.iloc[t,0] = time.time() - start_time
    #print("--- %s seconds ---" % (t))
    #record_ngmm.iloc[t,2] = ngmm.bic(data)
    record_ngmm.iloc[t,1] = ngmm.score(data)
    record_ngmm.iloc[t,2] = final_k#trav.n_components
    #print(ngmm_ll)
    
    
    print('naive GMM with fix k.')
    start_time = time.time()
    ngmm = GaussianMixture(trav.n_components,covariance_type='diag').fit(data)
    record_ngmm_fix_k.iloc[t,0] = time.time() - start_time
    #print("--- %s seconds ---" % (t))
    #record_ngmm.iloc[t,2] = ngmm.bic(data)
    record_ngmm_fix_k.iloc[t,1] = ngmm.score(data)
    record_ngmm_fix_k.iloc[t,2] = trav.n_components
    #print(ngmm_ll)
    
    
    print('dpgmm.')
    start_time = time.time()
    dpgmm = BayesianGaussianMixture(n_components=max_cluster_num,max_iter=500).fit(data)
    record_dpgmm.iloc[t,0] = time.time() - start_time
    record_dpgmm.iloc[t,1] = dpgmm.score(data)
    record_dpgmm.iloc[t,2] = len(dpgmm.weights_)
    
    
db_summary = pd.concat([record_sort,record_gmm,record_gmm_fix_k,record_ngmm,record_ngmm_fix_k,record_dpgmm])
db_summary['DB'] = name
db_summary['method'] = ['CITE-sort']*record_sort.shape[0] + ['GMM']*record_gmm.shape[0] + ['GMM_fixk']*record_gmm_fix_k.shape[0] + \
['nGMM']*record_ngmm.shape[0] + ['nGMM_fixk']*record_ngmm_fix_k.shape[0] + ['dpgmm']*record_dpgmm.shape[0]

db_summary.to_csv(savepath+'/record_'+name+'.csv')




'''
record_full[name] = db_summary


record_full_alldb = pd.concat([record_full[name] for name in namelist])
record_full_alldb.to_csv(savepath+'/record_8DBs.csv')



plt.figure(figsize=(6,3), dpi=96)
ax = sns.barplot(x='DB', y='time', hue='method', data=record_full_alldb)
plt.ylabel('Time (s)')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(savepath+'./time.png')



plt.figure(figsize=(6,3), dpi=96)
ax = sns.barplot(x='DB', y='ll', hue='method', data=record_full_alldb)
plt.ylabel('log-likelihood')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(savepath+'./ll.png')
'''