# T-cell-effectorness
This repository contains all the codes used for the downstream data analysis steps in the paper "Single-cell transcriptomics identifies an effectorness gradient shaping the response of CD4+ T cells to cytokines"
For a detailed explanation of the aims and theoretical bases of these codes, please refer to our study in Nature Communications:

https://www.nature.com/articles/s41467-020-15543-y

All codes are written in R (either as R files or R markdowns) and have been organised by data type (see below for a detailed description of the data analysis steps contained in each directory).
For further information about these codes, please either raise an issue or email ecg@sanger.ac.uk

## Bulk RNA-sequencing
Contains all the codes used for analysing bulk RNA expression in T cells upon stimualtion in the presence of cytokines. The analysis steps contained in the repository include:

1) Integration of the counts table and metadata information into a single SummarizedExperiment object
2) Exploratory data analysis (i.e. normalisation and log-transformation, batch regression, pca and visualisation)
3) Differential expression analysis between resting and activated T cells, as well as between T cells activated in the presence or absence of cytokines

## Proteomics
Contains all the codes used for analysing bulk protein expression (i.e. LC-MS/MS data) in T cells upon stimulation in the presence of cytokines. The analysis steps contained in the repository include:

1) Integration of the relative protein abundances table and metadata information into a single SummarizedExperiment object
2) Exploratory data analysis (i.e. batch regression, pca and visualisation)
3) Differential expression analysis between resting and activated T cells, as well as between T cells activated in the presence or absence of cytokines

## Multiomic analysis of RNA and protein expression
Contains all the codes used for integrating RNA-seq and quantitative proteomics data. The analysis steps contained in this repository include:

1) Integration of RNA counts and relative protein abundances into a single "multiomics" expression marix
2) A multiomic differential expression analysis which takes into account both RNA and protein expression values for each gene. This analysis is based on the f-divergence cutoff index (fCI) method.
3) Identification of cell type specific "proteogenomic" signatures. That is, genes which are expressed specifically upon stimultion with a given cytokine, both at the RNA and protein level. The functions for this analysis are available in the R package 'proteogenomic' (https://github.com/eddiecg/proteogenomic)

## Single-cell RNA-sequencing
Contains all the codes used for analysing single-cell gene expression (i.e. 3' 10X data) in T cells upon stimulation in the presence or absence of cytokines. The analysis steps in this repository include:

1) Deconvolution of cells of different individuals based on natural genetic variation. This step of the analysis is based on the Cardelino (now Vireo) algorithm (https://github.com/single-cell-genetics/vireo). 
2) Integration of single-cell RNA counts and cell annotations from different sample into a single expression matrix
3) Exploratory data analysis of over 40,000 single T cells from different conditions (e.g. normalisation, log-transformation, dimensionality reduction and visualisation)
4) Unsupervised clustering of resting T cells
5) Unsupervised clustering of T cells stimulated in the presence or different cytokines and identification of cytokine-specific cell states
6) Comparison of T cell response to different cytokines at the level of cell populations. This analysis is based on the single-cell implementation of the UniFrac method (https://github.com/liuqivandy/scUnifrac)
7) Ordering of T cells into pseudotime trajectories, which reveals the existence of a "T cell effectorness" gradient
8) Analysis of the relationship between T cell effectorness and TCR clonality. This analysis is based on publicly available single-cell paired TCR-sequences
9) Comparison of T cell effectorness between resting and stimulated cells based on canonical correlation analysis
10) Modelling of changes in gene expression as a function of T cell effectorness and cytokines  
