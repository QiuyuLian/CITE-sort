# CITE-sort

A new computational framework for auto-gating and robust segregation of artificial cell types in CITE-seq.

## Description

CITE-sort conducts auto-gating with CITE-seq ADT data using recursive Gaussian Mixture Model. It is robust against artificial cell types that stem from multiplets. CITE-sort also provides concrete explanations of its internal decision process by constructing a biologically meaningful Single Cell Taxonomy tree.  

Below shows an example Single Cell Taxonomy tree constructed by CITE-sort from an in-house PBMC dataset. Each node represents a subpopulation. The title of each inner node represents the selected surface markers subspace. Red and blue colors represent the two component complexes for subdivision. Edges are colored according to their corresponding component complexes. Leaf nodes are hand-curated and are annotated with domain knowledge. Cell types that should not exist are labeled as suspect _artificial cell type_ (ACT) clusters. Suspect ACT clusters are characterized by their population percentages in the overall dataset (denoted by ‘prop’) and their multi-sample multiplets percentages (denoted by ‘MSM’). Abbreviations: iNK: intermediate NK cells; mNK: vast majority of NK cells; C-mono: classical monocytes; NC-mono: non-classical monocytes; mDC: myeloid DC; DNT: double negative T cells.

<img src="readme_figs/taxonomy.png" alt="taxonomy" style="zoom:67%;" />

## Usage

### Input

The input of CITE-sort should be a csv file with CLR normalized CITE-seq ADT data (row: droplet/sample, col: ADT/feature). 

### Run

`python CITEsort.py ADT_clr_file -c 0.1 -o ./CITEsort_out`

- -c, cutoff, the similarity threshold of merging Gaussian components; the default is 0.1. It should be a real value between 0 and 1. The bigger value leads to split more aggressively, and ends in a more complicated tree.
- -o, output, the path to save ouput files. If not specified, CITE-sort will create a folder "./CITEsort_out" in the current directory.

### Outputs

- tree.pdf, the vasualized Single Cell Taxonomy tree of input dataset created by CITE-sort.
- leaf_labels.csv, the labels of each droplets in Single Cell Taxonomy tree.
- tree.pickle, the tree structure recording the main clusteirng infromation of input dataset.
- tree.dot, the auxiliary file to plot the tree.

## Examples

We provide 4 public CITE-seq datasets as examples in "./data_10x":

- [10k PBMCs from a healthy donor](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3) 
- [10k cells from a MALT tumor](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/malt_10k_protein_v3)
- [5k PBMCs from a healthy donor (next GEM)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3_nextgem)
- [5k PBMCs from a healthy donor (v3 chemistry)](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3)

### Example Commond

The 10k PBMCs dataset is used as example. 

`python CITEsort.py ./data_10x/ADT_10k_PBMC_normalized.csv -o ./10x_10k_PBMC `

## Authors

Qiuyu Lian\*, Hongyi Xin\*, Jianzhu Ma, Yanda Li, Liza Konnikova, Wei Chen\#, Jin Gu\#,Kong Chen\#

## Maintainer

Qiuyu Lian, Hongyi Xin.



