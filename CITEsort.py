#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 23:49:04 2019

@author: lianqiuyu
"""

import pandas as pd
from Matryoshka import Matryoshka
from Visualize import visualize_tree
from BTreeTraversal import BTreeTraversal
import pickle
import argparse
import os

#from sys import argv

def main():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('data_path', help = "The input path of CLR normalized data in .csv files with row as sample, col as feature.")
    parser.add_argument('-c','--cutoff', help = "The cutoff for merging components (default 0.1). It shoube between 0 and 1. The bigger value leads to split more aggressively, and ends in a more complicated tree.")
    parser.add_argument("-o", "--output", help="The path for storing the tree structure and leaf labels. Requires a path argument.", type=str)

    args = parser.parse_args()
    data_path = args.data_path
    
    if args.output:
        output_path = args.output
    else:
        output_path = "./CITEsort_out"
    
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    if args.cutoff:
        merge_cutoff = float(args.cutoff)
    else:
        merge_cutoff = 0.1
    
    data = pd.read_csv(data_path,header = 0, index_col=0)
    tree = Matryoshka(data,merge_cutoff)
    visualize_tree(tree,data,output_path,'tree')
    
    f = open(output_path+'/tree.pickle','wb')
    pickle.dump(tree,f)
    f.close()
    
    traversal = BTreeTraversal(tree)
    leaves_labels = traversal.get_leaf_label()
    leaves_labels.to_csv(output_path + '/leaf_labels.csv')


if __name__ == "__main__":
    main()
    
    