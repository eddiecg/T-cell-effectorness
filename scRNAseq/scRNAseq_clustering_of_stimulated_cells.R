### UNSUPERVISED CLUSTERING OF STIMULATED CD4+ T CELLS ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code performs unsupervised clustering and cluster annotation of stimulated T cells.

# LOADING LIBRARIES
library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(reshape2)

# LOADING DATA
## Loading a Seurat object containing the expression matrix and annotations for all cells in the study
## Please refer to "scRNAseq_exploratory_data_analysis.R" to see how this object was generated
T.cells <- readRDS("allCells_Seurat.rds")

# SUBSETTING DATA
## Keeping only resting T cells (i.e. unstimulated or UNS)
stim.cells <- SubsetData(T.cells, cells.use=T.cells@cell.names[T.cells@meta.data$cytokine.condition != "UNS"],subset.raw = T)

# REDUCING DIMENSIONS
## Identifying highly variable genes
stim.cells <- ScaleData(stim.cells, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))
stim.cells <- FindVariableGenes(stim.cells,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)

## Performing dimensionality reduction and embedding
stim.cells <- RunPCA(stim.cells, pc.genes = stim.cells@var.genes[!stim.cells@var.genes %in% cc.genes], do.print = FALSE, pcs.compute = 50)
PCElbowPlot(stim.cells, num.pc = 50)

stim.cells <- RunUMAP(stim.cells, dims.use=1:30, reduction.use = "pca")

# CLUSTERING DATA
## Clustering cells with the Louvain algorithm
stim.cells <- FindClusters(stim.cells, reduction.type = "pca", dims.use = 1:30, resolution = 1)

## Identifying markers of each cluster
stim.cells <- FindAllMarkers(stim.cells, only.pos=FALSE, min.pct=0.1, logfc.threshold=0.25)

## Annotating clusters based on their gene markers
### Adding a 'cluster_id' variable to the metadata and setting it to a dummy value
stim.cells@meta.data$cluster_id <- "Unidentified"

### Adding the correct 'cluster_id' annotation for each cell based on the identified marker genes, the clustering results and the respective cell type and cytokine condition annotations
stim.cells@meta.data[stim.cells@meta.data$res.1==0,]$cluster_id <- "TN (Th17)"
stim.cells@meta.data[stim.cells@meta.data$res.1==1,]$cluster_id <- "TN (iTreg)"
stim.cells@meta.data[stim.cells@meta.data$res.1==2,]$cluster_id <- "TCM2 (Th17/iTreg)"
stim.cells@meta.data[stim.cells@meta.data$res.1==3,]$cluster_id <- "TN (Th2)"
stim.cells@meta.data[stim.cells@meta.data$res.1==4,]$cluster_id <- "TN (Th17/iTreg)"
stim.cells@meta.data[stim.cells@meta.data$res.1==5,]$cluster_id <- "TCM2 (Th0)"
stim.cells@meta.data[stim.cells@meta.data$res.1==6 & stim.cells@meta.data$cell_type=="Memory" & !(stim.cells@meta.data$cytokine.condition %in% c("Th17","iTreg")),]$cluster_id <- "TEM (Th0)"
stim.cells@meta.data[stim.cells@meta.data$res.1==6 & stim.cells@meta.data$cell_type=="Memory" & (stim.cells@meta.data$cytokine.condition %in% c("Th17","iTreg")),]$cluster_id <- "TEM (Th17/iTreg)"
stim.cells@meta.data[stim.cells@meta.data$res.1==6 & stim.cells@meta.data$cell_type!="Memory",]$cluster_id <- "TN (Th0)"
stim.cells@meta.data[stim.cells@meta.data$res.1==7,]$cluster_id <- "TEM (Th17/iTreg)"
stim.cells@meta.data[stim.cells@meta.data$res.1==8,]$cluster_id <- "TCM1 (Th0)"
stim.cells@meta.data[stim.cells@meta.data$res.1==9,]$cluster_id <- "TN (Th0)"
stim.cells@meta.data[stim.cells@meta.data$res.1==10,]$cluster_id <- "Mitotic"
stim.cells@meta.data[stim.cells@meta.data$res.1==11 & (stim.cells@meta.data$cytokine.condition=="Th17" | stim.cells@meta.data$cytokine.condition=="iTreg"),]$cluster_id <- "TEMRA (Th17/iTreg)"
stim.cells@meta.data[stim.cells@meta.data$res.1==11 & !(stim.cells@meta.data$cytokine.condition=="Th17" | stim.cells@meta.data$cytokine.condition=="iTreg"),]$cluster_id <- "TEMRA (Th0)"
stim.cells@meta.data[stim.cells@meta.data$res.1==12,]$cluster_id <- "IFN.high"
stim.cells@meta.data[stim.cells@meta.data$res.1==13,]$cluster_id <- "HSP.high"
stim.cells@meta.data[stim.cells@meta.data$res.1==14,]$cluster_id <- "TCM1 (Th17/iTreg)"
stim.cells@meta.data[stim.cells@meta.data$res.1==15,]$cluster_id <- "nTreg (Th0)"
