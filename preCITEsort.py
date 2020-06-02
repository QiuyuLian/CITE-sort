#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 19:25:32 2020

@author: lianqiuyu
"""

import pandas as pd
import argparse
import os
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import sys


parser = argparse.ArgumentParser()
parser.add_argument('data_path', help = "The input path of CLR normalized data in .csv files with row as sample, col as feature.")
parser.add_argument('-o', '--output', type=str, default='./CITEsort_out',help='Path to save output files.')
parser.add_argument('--CLR', action='store_true', default=False, help='Input is raw counts. Transform counts into CLR format.')

args = parser.parse_args()
data_path = args.data_path

if not os.path.exists(data_path):
    print('Error: input file does not exist. Please check.')
    sys.exit(0)
    
if args.output:
    output_path = args.output
else:
    output_path = "./CITEsort_out"

if not os.path.exists(output_path):
    os.mkdir(output_path)

print('read data.')
data = pd.read_csv(data_path,header=0,index_col=0)
dataplot = data

if args.CLRTransder:
    print('perform CLR transformation on raw counts.')
    data_clr = np.apply_along_axis(lambda x: np.log(x+1) - np.mean(np.log(x+1)),0,data)
    data_clr = pd.DataFrame(data_clr,index=data.index,columns = data.columns)
    data_clr.to_csv(output_path+'/data_clr.csv')
    dataplot = data_clr

print('plot histgrams of all markers in CLR format.')
plt.figure(figsize=(12,2*np.ceil(data.shape[1] / 5)), dpi=96)
plt.style.use('seaborn-white')
for i in range(dataplot.shape[1]):
    ax = plt.subplot(int(np.ceil(dataplot.shape[1] / 5)),5,i+1)
    sns.distplot(dataplot.iloc[:,i].values,kde_kws={'bw':0.2})
    plt.yticks([0,1])
    plt.title(dataplot.columns[i],fontsize=15)
    if i%5 == 0:
        plt.ylabel('Density',fontsize=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')

plt.suptitle('DB: '+str(dataplot.shape[1])+' ADTs,'+str(dataplot.shape[0])+' droplets',fontsize=15)    
plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9, hspace=0.6,wspace=0.15)
#plt.subplots_adjust(top=0.85)
#plt.savefig('./PBMC_16k/marker_hist.png')
plt.savefig(output_path+'/data_cls_hist.png')
plt.clf()
#plt.show()

    
    