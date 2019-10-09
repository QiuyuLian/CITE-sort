#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:51:16 2019

@author: lqyair
"""

#import pandas as pd
import numpy as np
#from BTreeTraversal import BTreeTraversal
from matplotlib import pyplot as plt
from scipy import stats
import pandas as pd

#node = traversal.get_node(0)
#nodename = traversal.nodename[0]

def visualize_node(data,node,nodename,**plot_para):
    
    # plot_para: savefig, outpath, 
    savefig = plot_para.get('savefig',False)
    savepath = plot_para.get('savepath','.')
    savename = plot_para.get('savename','.')
    
    current_indices = node.indices
    node_data = data.loc[current_indices,:]
    
    plt.figure(figsize=(12,((data.shape[1]-1)//5+1)*2), dpi=96)
    plt.style.use('seaborn-white')
    #ax.tick_params(axis='both', which='major', labelsize=10)

    
    if node.key == ('leaf',) and node_data.shape[0] <= 20 :
        markers = node_data.columns.values.tolist()
        for i in range(len(markers)):
            X = node_data.loc[:,markers[i]].values.reshape(-1, 1)
            plt.subplot( (len(markers)-1)//5+1,5,i+1 )
            plt.hist(X,bins=30, density = True, color = "lightblue")
            plt.ylabel('density',fontsize=10)
            plt.title( markers[i],fontsize=12)

    else:
        all_clustering = node.all_clustering_dic[1]
        markers = list(all_clustering.keys())
        
        for i in range(len(markers)):
            
            X = node_data.loc[:,markers[i]].values.reshape(-1, 1)
            
            plt.subplot( (len(markers)-1)//5+1,5,i+1 )

            bins = np.linspace(min(X),max(X),500)
            cols = ['r','g','b','c','m','y','darkorange','lightgreen','lightpink','darkgray']
    
            bp_ncluster = int(all_clustering[markers[i]]['bp_ncluster'])
            mp_ncluster = 1 # default
            weights = all_clustering[markers[i]]['bp_pro']
            means = all_clustering[markers[i]]['bp_mean']
            sigmas = np.sqrt(all_clustering[markers[i]]['bp_Sigma'])
            
            y = np.zeros((len(bins),bp_ncluster))
            
            for k in range(bp_ncluster):
                y[:,k] = (weights[k] * stats.norm.pdf(bins, means[k], sigmas[k]))[:,0]
                plt.plot(bins,y[:,k],linewidth=0.6,color='black')

            if bp_ncluster > 1:
                mp_ncluster = all_clustering[markers[i]]['mp_ncluster']
                mergedtonumbers = all_clustering[markers[i]]['mergedtonumbers']
                
                for k in range(mp_ncluster):
            
                    merged_idx = [idx for idx,val in enumerate(mergedtonumbers) if val == k]
                    y_merged = np.apply_along_axis(sum,1,y[:,merged_idx])
        
                    plt.plot(bins,y_merged,cols[k],linewidth=2,linestyle='-.')
                    
            subfig_title = '_'.join(markers[i])+' ('+str(mp_ncluster)+'|'+str(bp_ncluster)+') ' + str(round(all_clustering[markers[i]]['similarity_stopped'],2))
            
            if markers[i] == node.key:
                plt.title( subfig_title,fontsize=12,color='red')
            else: 
                plt.title( subfig_title,fontsize=12,color='darkgrey' if mp_ncluster <= 1 else 'black')
                
            plt.hist(X,bins=30, density = True, color = "lightblue")
            plt.ylabel('density',fontsize=10)
    
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9, hspace=0.4,wspace=0.45)
    plt.suptitle(nodename+' | '+str(len(current_indices))+' cells',fontsize=15,color="darkblue")
    plt.subplots_adjust(top=0.8)
    #plt.savefig(savepath+'/visualize_node.png')
    if savefig == True:
        plt.savefig(savepath+'/'+savename+'_'+nodename+'.png') 
    plt.show()   
    
    




#import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

def visualize_pair(data,node,nodename,**plot_para):
    
    # plot_para: savefig, outpath, 
    savefig = plot_para.get('savefig',False)
    savepath = plot_para.get('savepath','.')
    savename = plot_para.get('savename','.')
    
    all_clustering = node.all_clustering_dic[2]
    marker_pairs = list(all_clustering.keys())
    current_indices = node.indices   

    plt.figure(figsize=(12,((len(marker_pairs)-1)//5+1)*2.5), dpi=96)
    sns.set_style("white")
    
    for i in range(len(marker_pairs)):
    
        marker1,marker2 = marker_pairs[i]
        X1 = data.loc[current_indices, marker1]
        X2 = data.loc[current_indices, marker2]
        
        bp_clustering = all_clustering[marker_pairs[i]]['bp_clustering']
        mp_clustering = all_clustering[marker_pairs[i]]['mp_clustering']
        
        mp_ncluster = all_clustering[marker_pairs[i]]['mp_ncluster']
        bp_ncluster = all_clustering[marker_pairs[i]]['bp_ncluster']
    
        data_pair = pd.DataFrame({marker1:X1,marker2:X2,
                              'bp':bp_clustering,
                              'mp':mp_clustering},index=node.indices)

        plt.subplot( (len(marker_pairs)-1)//5+1,5,i+1 )
        
        #shapes = ['s','X','+']
        #markers = dict(zip(np.unique(mp_clustering),[shapes[idx] for idx in range(mp_ncluster)]))
        sns.scatterplot(x=marker1, y=marker2,hue="bp",style="mp",
                        data=data_pair,s=15,legend=False);

        marker_pair_joint = marker_pairs[i][0]+'_'+marker_pairs[i][1]
        subfig_title = marker_pair_joint+' ('+str(mp_ncluster)+'|'+str(bp_ncluster)+') ' + str(round(all_clustering[marker_pairs[i]]['similarity_stopped'],2))
        
        if marker_pairs[i] == node.key:
            plt.title( subfig_title,fontsize=12,color='red')
        else: 
            plt.title( subfig_title,fontsize=12,color='darkgrey' if mp_ncluster <= 1 else 'black')
            
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9, hspace=0.6,wspace=0.45)
    plt.suptitle(nodename+' | '+str(len(current_indices))+' cells',fontsize=15,color="darkblue")
    plt.subplots_adjust(top=0.85)
    #plt.savefig(savepath+'/visualize_node.png')
    if savefig == True:
        plt.savefig(savepath+'/'+savename+'_'+nodename+'.png') 
    plt.show() 
    
  




    
from subprocess import call
#from IPython.display import Image
#import pandas as pd
import matplotlib
#import numpy as np

def visualize_tree(root,data,outpath,filename):
    """write tree structure into .dot and .png files."""
    
    # open a file, and design general format
    tree_dot = open(outpath+'/'+filename+'.dot','w') 
    tree_dot.writelines('digraph Tree {')
    tree_dot.writelines('node [shape=box, style="filled, rounded", color="black", fontname=helvetica] ;')
    tree_dot.writelines('edge [fontname=helvetica] ;')


    #tree_dot = _write_tree_bfs(root,tree_dot)
        # Base Case 
    if root is None: 
        return
    
    
    # Create an empty queue for level order traversal 
    queue = [] 
    nodelist = []
    idxStack = []
    
    tot_cells = len(root.indices)
    #means_in_root = root.marker_summary['mean']
    #stds_in_root = root.marker_summary['std']
    means_in_root = data.mean(axis = 0) 
    stds_in_root = data.std(axis = 0)
    markers = means_in_root.index.values.tolist()
    
    # auxiliary parameters for color display
    branch_col = pd.Series({1:'#ffccccff',2:'#ffff99ff',3:'#CC99CC',4:'#99CCFF'})   
    leaf_col = matplotlib.colors.Normalize(vmin=0, vmax=np.log(tot_cells))
    
    node = root
    
    # Enqueue Root and initialize height 
    queue.append(node) 
    
    i = 0
    #print(str(i)+'_'+root.key)
    all_clustering = node.all_clustering_dic[len(node.key)]
    bp_ncluster = all_clustering[node.key]['bp_ncluster']
    mp_ncluster = all_clustering[node.key]['mp_ncluster']
    tree_dot.writelines(str(i)+' [label="'+str(i)+'_'+'_'.join(node.key)+ \
                        '\\nNum: '+str(len(node.indices))+ \
                        '\\n('+str(mp_ncluster)+'|'+str(bp_ncluster)+')",fillcolor="#ff9966ff",fontsize=25];')  
    nodelist.append(node.key)
    idxStack.append(i)
    
    while(len(queue) > 0): 
        # Print front of queue and remove it from queue 
        node = queue.pop(0) 
        idx = idxStack.pop(0)
                
        # left child 
        if node.left is not None: 
            nodelist.append(node.left.key)
            queue.append(node.left)
            i = i + 1
            idxStack.append(i)
            #print(str(i)+'_'+node.left.key)
            
            percent = str(round(len(node.left.indices)/tot_cells*100,2))+'%'
            mean_temp = data.loc[node.left.indices,:].mean(0) 
            
            if node.left.key == ('leaf',):
                # left leaf node        
                temp = (mean_temp - means_in_root)/stds_in_root
                offset_in_leaf = '\n' + markers[0]+': '+str(round(temp[markers[0]],2))+')'
                for k in range(1,len(markers)):
                    offset_in_leaf = offset_in_leaf + '\n' + markers[k]+': '+ str(round(temp[markers[k]],2))
                
                col =  matplotlib.colors.to_hex(matplotlib.cm.Greens(leaf_col(np.log(len(node.left.indices)))))
                tree_dot.writelines(str(i)+' [label="'+str(i)+'_'+'_'.join(node.left.key)+'\\n'+ \
                                    str(len(node.left.indices))+ ' ('+percent+')\\n'+ \
                                    offset_in_leaf+'",fillcolor="'+col+'",fontsize=20];')
            else:
                # left branch node
                all_clustering = node.left.all_clustering_dic[len(node.left.key)]
                bp_ncluster = all_clustering[node.left.key]['bp_ncluster']
                mp_ncluster = all_clustering[node.left.key]['mp_ncluster']
                
                tree_dot.writelines(str(i)+' [label="'+str(i)+'_'+'_'.join(node.left.key)+'\\n'+ \
                                    str(len(node.left.indices))+' ('+percent+')\\n'+ \
                                    '('+str(mp_ncluster)+'|'+str(bp_ncluster)+')",fillcolor="'+branch_col[len(node.left.key)]+'",fontsize=25];')

            # edge from parent to left node
            offset = ''
            for m in nodelist[idx]:
                val = (mean_temp[m] - means_in_root[m])/stds_in_root[m]
                offset = offset + str(round(val,2))+'\n'
            #print(str(idx)+'->'+str(i))
            tree_dot.writelines(str(idx)+' -> '+str(i)+ ' [labeldistance=3, label = "'+offset+'",fontsize=25, color='+['black','red'][node.where_dominant=='left']+\
                                ', style='+['solid','bold'][node.where_dominant=='left']+'];')

        # right child 
        if node.right is not None: 
            nodelist.append(node.right.key)
            queue.append(node.right) 
            i = i + 1
            idxStack.append(i)
            #print(str(i)+'_'+node.right.key)
            
            percent = str(round(len(node.right.indices)/tot_cells*100,2))+'%'
            mean_temp = data.loc[node.right.indices,:].mean(0) 

            if node.right.key == ('leaf',):
                # right leaf node
                temp = (mean_temp - means_in_root)/stds_in_root
                offset_in_leaf = '\n' + markers[0]+': '+str(round(temp[markers[0]],2))
                for k in range(1,len(markers)):
                    offset_in_leaf = offset_in_leaf + '\n' + markers[k]+': '+ str(round(temp[markers[k]],2))

                col =  matplotlib.colors.to_hex(matplotlib.cm.Greens(leaf_col(np.log(len(node.right.indices)))))
                tree_dot.writelines(str(i)+' [label="'+str(i)+'_'+'_'.join(node.right.key)+'\\n'+ \
                                    str(len(node.right.indices))+ ' ('+percent+')'+'\\n'+ \
                                    offset_in_leaf+'",fillcolor="'+col+'",fontsize=20];')

            else:
                # right branch node
                all_clustering = node.right.all_clustering_dic[len(node.right.key)]
                bp_ncluster = all_clustering[node.right.key]['bp_ncluster']
                mp_ncluster = all_clustering[node.right.key]['mp_ncluster']
                
                tree_dot.writelines(str(i)+' [label="'+str(i)+'_'+'_'.join(node.right.key)+'\\n'+ \
                                    str(len(node.right.indices))+' ('+percent+')\\n'+ \
                                    '('+str(mp_ncluster)+'|'+str(bp_ncluster)+')",fillcolor="'+branch_col[len(node.right.key)]+'",fontsize=25];')

            # edge from parent to right node
            offset = ''
            for m in nodelist[idx]:
                val = (mean_temp[m] - means_in_root[m])/stds_in_root[m]
                offset = offset + str(round(val,2))+'\n'
            #print(str(idx)+'->'+str(i))
            tree_dot.writelines(str(idx)+' -> '+str(i)+' [labeldistance=3, label = "'+offset+'",fontsize=25, color='+['black','red'][node.where_dominant=='right']+ \
                                ', style='+['solid','bold'][node.where_dominant=='right']+'];')
    
    # main body is completed
  
    tree_dot.writelines('}')
    tree_dot.close()

    # Convert to png using system command (requires Graphviz)
    call(['dot', '-Tpdf', outpath+'/'+filename+'.dot', '-o', outpath+'/'+filename+'.pdf', '-Gdpi=100'])
    
    # Display in jupyter notebook
    #Image(filename = outpath+'/GatingTree.png')


