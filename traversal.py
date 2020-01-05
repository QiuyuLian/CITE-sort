#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 20:58:29 2019

@author: lianqiuyu
"""

import pandas as pd
import numpy as np
#from matplotlib import pyplot as plt
#from Visualize import visualize_node,visualize_pair

class Traversal:
    
    def __init__(self,tree,method='bfs',nodelist=None,nodename=None,posterior=None,leaf_label=None,ll_tot =None,prop=None):
        
        #print('initializing...')
        
        self.tree = tree
        self.method = method
        if self.method == 'bfs':
            self.nodelist = self.levelOrderTraversal()
        if self.method == 'dfs':
            self.nodelist = self.preorderTraversal()
        
        nodename_temp = ['_'.join(x.key) for x in self.nodelist]
        self.nodename = [str(i)+'_'+nodename_temp[i] for i in range(len(nodename_temp))]
        self.posterior, self.leaf_label,self.ll_tot, self.prop = self.summarize()
         
        
      
    def summarize(self):
        #print('summarizing...')
        ll_tot = 0
        leafnodename = [x for x in self.nodename if x.split('_')[1] == 'leaf']
        posterior = pd.DataFrame(np.zeros([len(self.nodelist[0].sample_weights),len(leafnodename)]),columns=leafnodename,index=self.nodelist[0].sample_weights.index)
        prop = pd.Series(np.zeros(len(leafnodename)),index=leafnodename)
        
        for leaf in leafnodename :
            nidx = self.nodename.index(leaf)
            posterior.loc[:,leaf] = self.nodelist[nidx].sample_weights
            ll_tot += self.nodelist[nidx].ll_tot
            prop[leaf] = self.nodelist[nidx].prop
            
        leaf_label = posterior.idxmax(1)

        return posterior, leaf_label, ll_tot, prop
        

    
    def get_node(self,nodeID):
        return self.nodelist[nodeID]
 



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
















