# CITE-sort

CITE-sort is an interpretable clustering framework for CITE-seq datasets based on recursive Gaussian Mixture Model. 

## Description

CITE-sort conducts auto-gating with CITE-seq ADT data. It is robust against artificial cell types that stem from multiplets. CITE-sort also generates biologically meaningful interpretations to its clustering results by constructing a Single Cell Taxonomy tree (Below shows an eaxmple result).

<img src="/Users/qiuyulian/Documents/project/Matryoshka/readme_figs/taxonomy.png" alt="taxonomy" style="zoom:60%;" />

We call cell types that truly exist as *biological cell types* (BCT) and cell types created by CITE-seq multiplets as *artificial cell types* (ACT).   ACTs induce high imbalance of both cluster sizes and clustering coefficients. It seriously influence the performance and efficienty of in clustering. 

![imbalance](/Users/qiuyulian/Documents/project/Matryoshka/readme_figs/ACTimbalance.png)

CITE-sort addresses this problem by transforming the original high-dimensional, highly imbalanced, many-class clustering problem into solving multiple low-dimensional, fewer-class clustering sub-problems.   

## Usage

### Input

The input of CITE-sort 

### Run

`python CITEsort.py <data_path>`

### Outpus





