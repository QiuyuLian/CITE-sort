#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 23:44:26 2020

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
'''
namelist = ['PBMC_1k','PBMC_1k_b','PBMC_2k', 'PBMC_5k', 'PBMC_8k', 'MALT_8k', 'CBMC_8k','PBMC_16k']
datapath = './datasets/data_for_performance'
savepath = './performance'

record_full = {}
for name in namelist:
    db_summary = pd.read_csv(savepath+'/record_'+name+'.csv',header=0,index_col=0)
    record_full[name] = db_summary


record_full_alldb = pd.concat([record_full[name] for name in namelist])
record_full_alldb.to_csv(savepath+'/record_8DBs.csv')

'''
record_full_alldb = pd.read_csv('./performance/record_8DBs.csv',header=0,index_col=0)

temp = record_full_alldb.loc[record_full_alldb['method']!='GMM_fixk',]
record_plot = temp.loc[temp['method']!='nGMM_fixk',:]

record_plot['time'] = record_plot['time']/60





plt.figure(figsize=(8,3), dpi=96)
ax = sns.barplot(x='DB', y='time', hue='method', data=record_plot)
plt.ylabel('Time (min)',fontsize=15)
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('./performance/time.pdf')
plt.show()


record_plot['ll'] = - record_plot['ll']

plt.figure(figsize=(8,3), dpi=96)
ax = sns.barplot(x='DB', y='ll', hue='method', data=record_plot)
plt.ylabel(' - log-likelihood',fontsize=15)
plt.xlabel('')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.savefig('./performance/ll.pdf')
plt.show()


