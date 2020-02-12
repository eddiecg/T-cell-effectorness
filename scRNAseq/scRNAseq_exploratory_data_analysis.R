### EXPLORATORY ANALYSIS OF SINGLE-CELL RNA-SEQ DATA ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code performs normalization, regression of covariates and scaling of scRNA-seq data. Next, it performs dimensionality reduction and visualisation.

# LOADING LIBRARIES
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)

# LOADING DATA
## Loading expression matrix from a directory containing the matrix, genes and barcodes file
#cells <- Read10X("./scRNAseq_data/")
expressionMatrix <- Read10X("~/cytoimmgen/phase_1/data/")

## Loading cell annotations
#metadata <- read.table("./scRNAseq_data/NCOMMS-19-7936188_metadata.txt", header=T, sep="\t")
metadata <- read.table("~/cytoimmgen/phase_1/data/NCOMMS-19-7936188_metadata.txt", header=T, sep="\t")

# CREATING SEURAT OBJECT
cells <- CreateSeuratObject(expressionMatrix, min.cells = 3, min.genes = 200, project = "Cytoimmgen_10X_optimization")

## Adding cell annotations as metadata
cells@meta.data <- metadata

# NORMALISING GENE EXPRESSION
cells <- NormalizeData(cells, normalization.method = "LogNormalize")

# REGRESSING COVARAITES
cells <- ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# FINDING HIGHLY VARIABLE GENES
cells <- FindVariableGenes(cells,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)

# REDUCING DIMENSIONS
cells <- RunPCA(cells, pc.genes = cells@var.genes, do.print = FALSE, pcs.compute = 50)

# UMAP
cells <- RunUMAP(cells, dims.use=1:30, reduction.use = "pca")

# VISUALISING DATA
DimPlot(cells, reduction.use = "umap", group.by = "cell.type", pt.size = 0.5)
DimPlot(cells, reduction.use = "umap", group.by = "donor.id", pt.size = 0.5)
DimPlot(cells, reduction.use = "umap", group.by = "Phase", pt.size = 0.5)
DimPlot(cells, reduction.use = "umap", group.by = "cytokine.condition", pt.size = 0.5)

# SAVING AS .RDS FILE
saveRDS(cells, "allCells_Seurat.rds")
