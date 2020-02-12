### UNSUPERVISED CLUSTERING OF RESTING CD4+ T CELLS ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code includes an exploratory analysis of scRNA-seq data of resting CD4+ T cells. Next, it performs unsupervised clustering and cluster annotation.

# LOADING LIBRARIES
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)

# LOADING DATA
## Loading a Seurat object containing the expression matrix and annotations for all cells in the study
## Please refer to "scRNAseq_exploratory_data_analysis.R" to see how this object was generated
T.cells <- readRDS("allCells_Seurat.rds")

# SUBSETTING DATA
## Keeping only resting T cells (i.e. unstimulated or UNS)
uns.cells <- SubsetData(T.cells, cells.use=T.cells@cell.names[T.cells@meta.data$cytokine.condition == "UNS"],subset.raw = T)

# EXPLORATORY ANALYSIS OF RESTING T CELLS
## Identifying highly variable genes
uns.cells <- ScaleData(uns.cells, vars.to.regress = c("donor.id","percent.mito","nUMI"))
uns.cells <- FindVariableGenes(uns.cells,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)

## Performing dimensionality reduction and embedding
uns.cells <- RunPCA(uns.cells, pc.genes = uns.cells@var.genes, do.print = FALSE, pcs.compute = 30)
PCElbowPlot(uns.cells, num.pc = 50)

uns.cells <- RunUMAP(uns.cells, dims.use=1:15, reduction.use = "pca")

# UNSUPERVISED CLUSTERING OF RESTING T CELLS
## Clustering cells with the Louvain algorithm
uns.cells <- FindClusters(object = uns.cells, reduction.type = "pca", dims.use = 1:15, resolution = 0.6, force.recalc=T,print.output=F)

## Identifying markers of each cluster
uns.markers <- FindAllMarkers(uns.cells, only.pos=FALSE, min.pct=0.1, logfc.threshold=0.25)

### Annotating clusters based on their gene markers
uns.cells@meta.data$cluster_id <- as.character(uns.cells@meta.data$cluster_id)

uns.cells@meta.data$cluster_id[uns.cells@meta.data$res.0.6==0] <- "TN (resting)"
uns.cells@meta.data$cluster_id[uns.cells@meta.data$res.0.6==1] <- "TCM (resting)"
uns.cells@meta.data$cluster_id[uns.cells@meta.data$res.0.6==2] <- "TEM (resting)"
uns.cells@meta.data$cluster_id[uns.cells@meta.data$res.0.6==3] <- "nTreg (resting)"
uns.cells@meta.data$cluster_id[uns.cells@meta.data$res.0.6==4] <- "TEMRA (resting)"

uns.cells@meta.data$cluster_id <- factor(uns.cells@meta.data$cluster_id)

### Setting the new cluster labels as cell identity variables
uns.cells@ident <- uns.cells@meta.data$cluster_id
names(uns.cells@ident) <- uns.cells@cell.names

# QUANTIFYING THE CONTRIBUTION OF DIFFERENT INDIVIDUALS TO EACH CLUSTER
## Estimating the proportion of cells in each cluster which belong to different individuals
props <- t(sapply(names(table(umap.res$donor.id)), FUN=function(donor){
  props <- table(umap.res[umap.res$donor.id==donor,]$cluster_id)
  props <- props/sum(props)
  return(props)
}))

props <- melt(props)
colnames(props) <- c("donor.id","cell_type","proportion")
props$cell_type <- factor(props$cell_type, levels=rev(c("T naive", "T CM", "T EM", "T EMRA", "nTreg")))

## Visualising the estimated proportions
ggplot(data=props, aes(x=donor.id)) + 
  geom_bar(aes(weight=proportion, fill=cell_type),position=position_stack()) + 
  xlab("Biological replicate") + 
  ylab("Proportion of cells") + 
  theme_test() + 
  coord_flip() + 
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12))
