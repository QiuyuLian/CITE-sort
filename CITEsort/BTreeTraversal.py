#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 20:26:57 2019

@author: lqyair
"""
import sys
sys.path.append("./CITEsort")

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from Visualize import visualize_node,visualize_pair

   
'''
from matplotlib import pyplot as plt
from matplotlib import cm,colors
from mpl_toolkits.axes_grid1 import axes_grid
#import seaborn as sns
import pdb
'''

class BTreeTraversal:
    
    def __init__(self,tree,method='bfs',nodelist=None,nodename=None,tree_summary=None,leaf_summary=None,ll=None,n_components=None):
        
        #print('initializing...')
        
        self.tree = tree
        self.method = method
        if self.method == 'bfs':
            self.nodelist = self.levelOrderTraversal()
        if self.method == 'dfs':
            self.nodelist = self.preorderTraversal()
        
        nodename_temp = ['_'.join(x.key) for x in self.nodelist]
        self.nodename = [str(i)+'_'+nodename_temp[i] for i in range(len(nodename_temp))]
        self.tree_summary, self.leaf_summary = self.summarize()
        if 'll' in self.tree.__dir__():
            self.ll = self.leaf_summary['ll'].sum()
            self.n_components = self.leaf_summary.shape[0]
        
      
    def summarize(self):
        if 'll' in self.tree.__dir__():
            tree_summary = pd.DataFrame({'Count':[len(x.indices) for x in self.nodelist],
                                         'Weight':[x.weight for x in self.nodelist],
                                         'Stop':[x.stop for x in self.nodelist],
                                         'll':[x.ll for x in self.nodelist]
                                         },index=self.nodename)
        else:
            tree_summary = pd.DataFrame({'Count':[len(x.indices) for x in self.nodelist] },index=self.nodename)  

        leaf_summary = tree_summary.loc[[x for x in self.nodename if x.split('_')[1]=='leaf'],:]
        leaf_summary = leaf_summary.sort_values(by='Count',ascending=False)
        
        return tree_summary,leaf_summary
    
    
    
    def get_ll(self):
        
        self.ll_tot = sum([x.ll for idx,x in enumerate(self.nodelist) if self.nodename[idx].split('_')[1]=='leaf'])
    
    
    def get_node(self,nodeID):
        return self.nodelist[nodeID]
 
    
    
    def generateLabel(self):
        """generate label file (binary matrix: Num.cell x Num.node, X_ij = 1 means cell i is attached to node j.)"""
        
        label = pd.DataFrame(np.zeros([len(self.tree.indices),len(self.nodename)]),index=self.tree.indices,columns=self.nodename)

        for i in range(len(self.nodename)):
            label.loc[self.nodelist[i].indices,self.nodename[i]] = 1

        return label
    
    
    
    def get_leaf_label(self):
        """generate label (one column, indicating which leaf cells are assigned.)"""
        label = pd.DataFrame({'GEM':self.tree.indices,'Label':[None]*len(self.tree.indices)},index=self.tree.indices)
        for i in range(len(self.nodename)):
            if self.nodename[i].split('_')[1] == 'leaf':
                label.loc[self.nodelist[i].indices,'Label'] = self.nodename[i]
                
        return label
    

    
    def plot_node(self,data, nodeID, viz_dim = 1, **plot_para):
        """plot the specified node (default: savefig=False,savepath='.')"""
        if viz_dim == 1:
            visualize_node(data, node = self.nodelist[nodeID], nodename = self.nodename[nodeID], **plot_para)
        if viz_dim == 2:
            visualize_pair(data, node = self.nodelist[nodeID], nodename = self.nodename[nodeID], **plot_para)

    
    
    def plot_leaf_size(self):
        
        leaf_size = self.leaf_summary['Count']
        leaf_prop = self.leaf_summary['Proportion']
        
        fig, ax1 = plt.subplots(1,2,figsize=(12,4))
        
        # plot number/proportion of cells in each leaf
        color = 'tab:red'
        ax1[0].set_xlabel('leaf',fontsize=20)
        ax1[0].set_ylabel('Proportion', color=color,fontsize=20)
        ax1[0].plot(range(len(leaf_prop)),leaf_prop, color=color,marker='o')
        if len(leaf_prop) >= 5:
            plt.xticks(np.arange(0, len(leaf_prop), len(leaf_prop)//5))
        else:
            plt.xticks(np.arange(0, len(leaf_prop), 1))
        ax1[0].tick_params(axis='y', labelcolor=color,labelsize=15)
        ax1[0].tick_params(axis='x', labelsize=15)
        ax1[0].set_title('Num. of cells in leaf',fontsize=20,pad=20)

        ax2 = ax1[0].twinx()  # instantiate a second axes that shares the same x-axis
        color = 'tab:blue'
        ax2.set_ylabel('Count', color=color,fontsize=20)  # we already handled the x-label with ax1
        ax2.plot(range(len(leaf_size)),leaf_size, color=color,marker='o')
        ax2.tick_params(axis='y', labelcolor=color,labelsize=15)

        
        # plot cumulative number/proportion of cells in each leaf
        color = 'tab:red'
        ax1[1].set_xlabel('leaf',fontsize=20)
        ax1[1].set_ylabel('Proportion', color=color,fontsize=20)
        ax1[1].plot(range(len(leaf_prop)),leaf_prop.cumsum(), color=color,marker='o')
        if len(leaf_prop) >= 5:
            plt.xticks(np.arange(0, len(leaf_prop), len(leaf_prop)//5))
        else:
            plt.xticks(np.arange(0, len(leaf_prop), 1))
        ax1[1].tick_params(axis='y', labelcolor=color,labelsize=15)
        ax1[1].tick_params(axis='x', labelsize=15)
        ax1[1].set_title('Cumulative num. of cell in leaf',fontsize=20,pad=20)

        ax2 = ax1[1].twinx()  # instantiate a second axes that shares the same x-axis
        color = 'tab:blue'
        ax2.set_ylabel('Count', color=color,fontsize=20)  # we already handled the x-label with ax1
        ax2.plot(range(len(leaf_size)),leaf_size.cumsum(), color=color,marker='o')
        ax2.tick_params(axis='y', labelcolor=color,labelsize=15)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()
        
    
    
    '''
    def track_marker(self,data,n_big_leaf,**plot_para):
        """track marker distributions in big leafs. (at most 12 leafs. default para: savefig=False,outpath='.')"""
        n_big_leaf = min(n_big_leaf,12) 
        
        savefig = plot_para.get('savefig',False)
        outpath = plot_para.get('outpath','.')
             
        big_leaf = self.leaf_summary.index.values.tolist()[0:n_big_leaf]
        markers = data.columns.values.tolist()
        node_plot = [self.nodename[0]] + big_leaf
        
        cmap = cm.get_cmap('Set3')
        col_dic = dict(zip(big_leaf,[colors.to_hex(cmap(i)) for i in range(len(big_leaf))]))
        col_dic[node_plot[0]] = '#999999' # col for all cells
        
        nrows = np.ceil(len(markers)/5)
        ncols = 5
        naxes = len(node_plot)
        f = plt.figure(figsize=(10, naxes))
        for i, m in enumerate(markers):
            ag = axes_grid.Grid(f, (nrows, ncols, i+1), (naxes, 1), axes_pad=0)
            for j in range(naxes):
                leaf_idx = int(node_plot[j].split('_')[0])
                ag[j].hist(data.loc[self.nodelist[leaf_idx].indices,m], 
                           color = col_dic[node_plot[j]], density = True, bins = 'auto')
                ag[j].axvline(0, linestyle='dashed', linewidth=2)
                ag[j].yaxis.set_ticks([])
                ag[j].xaxis.set_ticks([])
                if j%naxes == 0:
                    ag[j].set_title(markers[i],fontsize=15)
                if i%ncols == 0:
                    ag[j].set_ylabel(str(leaf_idx),fontsize=12)

        plt.subplots_adjust(wspace=0.1, hspace=0.3)
        if savefig == True:
            plt.savefig(outpath+'/track_marker_in_big_leafs.png')
        plt.show()
    '''
    
    
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






