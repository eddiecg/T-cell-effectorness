### ALIGNMENT OF RESTING AND TCR/CD28-STIMULATED CELLS ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code maps TCR/CD28-stimulated (i.e. Th0-stimulated) single-cells to their equivalents in the resting state
# This mapping is based on the alignment of both data sets into a common reduced space using cannonical correlation analysis (CCA)
# For more information on the details of CCA alignment, please refer to: https://doi.org/10.1038/nbt.4096

# LOADING LIBRARIES
library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(monocle)

# LOADING DATA
## Loading a Seurat object containing the expression matrix and annotations for all cells in the study
## Please refer to "scRNAseq_exploratory_data_analysis.R" to see how this object was generated
T.cells <- readRDS("allCells_Seurat.rds")

## Loading a list of cell cycle genes defined by the Regev lab
cc.genes <- readLines("regev_lab_cell_cycle_genes.txt")

# EXTRACTING RESTING AND TH0-STIMULATED CELLS
uns.cells <- SubsetData(T.cells, cells.use=T.cells@cell.names[T.cells@meta.data$cytokine.condition == "UNS"])
Th0 <- SubsetData(T.cells, cells.use=T.cells@cell.names[T.cells@meta.data$cytokine.condition == "Th0"])

# REMOVING CYCLING CELLS, AS WELL AS IFN AND HSP CLUSTERS
Th0 <- SubsetData(Th0, cells.use=Th0@cell.names[!Th0@meta.data$corrected_clusters %in% c("HSP.high","IFN.high", "Cycling.cells")])

## IDENTIFYING HIGHLY VARIABLE GENES
uns.cells <- FindVariableGenes(uns.cells,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)
Th0 <- FindVariableGenes(Th0,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)

# ALIGNING DATA SETS WITH CCA
## Identifying a union set of variable genes
union.genes <- setdiff(unique(c(uns.cells@var.genes, Th0@var.genes)),cc.genes)

## Performing CCA
cca.cells <- RunCCA(object=uns.cells, object2=Th0, genes.use=union.genes, num.cc=20)

## Discarding cells which remain poorly matched between data sets
cca.cells <- CalcVarExpRatio(object = cca.cells, reduction.type = "pca", grouping.var = "cytokine.condition", dims.use = 1:20)
cca.cells.all <- cca.cells
cca.cells <- SubsetData(object = cca.cells, subset.name = "var.ratio.pca", accept.low = 0.5)
rm(cca.cells.all)

## Aligning CCA subspaces
cca.cells <- AlignSubspace(object = cca.cells, reduction.type = "cca", grouping.var = "cytokine.condition", dims.align = 1:20)

## Running UMAP on the CCA dimensions
cca.cells <- RunUMAP(object = cca.cells, reduction.use = "cca.aligned", dims.use = 1:20)

# VISUALISING ALIGNED CELLS
## Extracting UMAP coordinates
umap.res <- data.frame(cbind(cca.cells@dr$umap@cell.embeddings, cca.cells@meta.data))

## Splitting coordiantes by cytokine condition
umap.stim <- umap.res[umap.res$cytokine.condition!="UNS",] 
umap.uns <- umap.res[umap.res$cytokine.condition=="UNS",]
umap.uns$cluster_id <- factor(umap.uns$cluster_id, levels=c("T naive","T CM","T EM","T EMRA","nTreg"))

## Visualising cells
cols <- c("#6baed6","#2171b5","#08519c","#08306b","#980043")
g1 <- ggplot(umap.stim, aes(x=UMAP1,y=UMAP2)) + geom_point(color="lightgrey", size=1, alpha=0.5) +
  stat_density2d(data=umap.uns[umap.uns$cluster_id=="T naive",], aes(x=UMAP1, y=UMAP2, alpha=..level..), fill=cols[4], size=2, bins=5, geom='polygon') +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  scale_alpha_continuous(range=c(0.2,0.8))
g2 <- ggplot(umap.stim, aes(x=UMAP1,y=UMAP2)) + geom_point(color="lightgrey", size=1, alpha=0.5) +
  stat_density2d(data=umap.uns[umap.uns$cluster_id=="T CM",], aes(x=UMAP1, y=UMAP2, alpha=..level..), fill=cols[4], size=2, bins=5, geom='polygon') +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  scale_alpha_continuous(range=c(0.2,0.8))
g3 <- ggplot(umap.stim, aes(x=UMAP1,y=UMAP2)) + geom_point(color="lightgrey", size=1, alpha=0.5) +
  stat_density2d(data=umap.uns[umap.uns$cluster_id=="T EM",], aes(x=UMAP1, y=UMAP2, alpha=..level..), fill=cols[4], size=2, bins=5, geom='polygon') +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  scale_alpha_continuous(range=c(0.2,0.8))
g4 <- ggplot(umap.stim, aes(x=UMAP1,y=UMAP2)) + geom_point(color="lightgrey", size=1, alpha=0.5) +
  stat_density2d(data=umap.uns[umap.uns$cluster_id=="T EMRA",], aes(x=UMAP1, y=UMAP2, alpha=..level..), fill=cols[4], size=2, bins=5, geom='polygon') +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  scale_alpha_continuous(range=c(0.2,0.8))
g5 <- ggplot(umap.stim, aes(x=UMAP1,y=UMAP2)) + geom_point(color="lightgrey", size=1, alpha=0.5) +
  stat_density2d(data=umap.uns[umap.uns$cluster_id=="nTreg",], aes(x=UMAP1, y=UMAP2, alpha=..level..), fill=cols[5], size=2, bins=5, geom='polygon') +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
  scale_alpha_continuous(range=c(0.2,0.8))
plot_grid(g1,g2,g3,g4, ncol=4,nrow=1)
