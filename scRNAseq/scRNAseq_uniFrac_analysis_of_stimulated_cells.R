### VISUALISATION AND SC-UNIFRAC ANALYSIS OF STIMULATED CD4+ T CELLS ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code includes an exploratory analysis of scRNA-seq data of stimulated CD4+ T cells. 
# Next, it performs dimensionality reduction and uses scUniFrac analysis to compare the location of cells stimulated with different cytokines in the reduced space.

# LOADING R LIBRARIES
library(Seurat)
library(scUnifrac)

# LOADING DATA
## Loading a Seurat object containing the expression matrix and annotations for all cells in the study
## Please refer to "scRNAseq_exploratory_data_analysis.R" to see how this object was generated
T.cells <- readRDS("allCells_Seurat.rds")

# SUBSETTING DATA
## Keeping only stimulated T cells
stim.cells <- SubsetData(T.cells, cells.use=T.cells@cell.names[T.cells@meta.data$cytokine.condition != "UNS"],subset.raw = T)

# EXPLORATORY DATA ANALYSIS OF STIMULATED T CELLS
## Identifying highly variable genes
stim.cells <- ScaleData(stim.cells, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))
stim.cells <- FindVariableGenes(stim.cells,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)

## Performing dimensionality reduction and embedding
stim.cells <- RunPCA(stim.cells, pc.genes = stim.cells@var.genes[!stim.cells@var.genes %in% cc.genes], do.print = FALSE, pcs.compute = 50)
PCElbowPlot(stim.cells, num.pc = 50)

stim.cells <- RunUMAP(stim.cells, dims.use=1:30, reduction.use = "pca")

# TESTING FOR DIFFERENCES IN "CELL POPULATION DIVERSITY" BETWEEN CONDITIONS
## Dividing cells into groups by cell type and cytokine condition
groups<-c(rep("Naive_Th0",sum(stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="Th0")),
          rep("Naive_Th2",sum(stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="Th2")),
          rep("Naive_Th17",sum(stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="Th17")),
          rep("Naive_iTreg",sum(stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="iTreg")),
          rep("Memory_Th0",sum(stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="Th0")),
          rep("Memory_Th2",sum(stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="Th2")),
          rep("Memory_Th17",sum(stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="Th17")),
          rep("Memory_iTreg",sum(stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="iTreg")))

## Loading cell cycle genes defined by the Regev lab
cc.genes <- readLines("regev_lab_cell_cycle_genes.txt")

## Computing pairwise UniFrac distances between all combinations of cell groups (i.e. distances between different cell type and cytokine condition combinations)
## For a more detailed explanation of scUniFrac, please refer to: https://doi.org/10.1371/journal.pbio.2006687
unifracs <- scUnifrac_multi(dataall=as.matrix(cbind(stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="Th0"],
                                                    stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="Th2"],
                                                    stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="Th17"],
                                                    stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Naive" & stim.cells@meta.data$cytokine.condition=="iTreg"],
                                                    stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="Th0"],
                                                    stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="Th2"],
                                                    stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="Th17"],
                                                    stim.cells@data[!rownames(stim.cells@data) %in% cc.genes,stim.cells@meta.data$cell_type=="Memory" & stim.cells@meta.data$cytokine.condition=="iTreg"])),
                            group=groups,
                            nDim = 4,
                            ncluster = 10,
                            normalize = F)
