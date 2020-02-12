### ORDERING OF T CELLS INTO EFFECTORNESS TRAJECTORIES ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code performs trajectory inference for groups of T cells exposed to different cytokines
# Using the monocle (v2) alogirthm, cells from each cytokine cytokine condition are ordered into a trajectory which starts with naive T cells and progresses towards highly effector memory T cells
# Each cell is assigned a pseudotime value based on its position along the trajectory
# Finally, pseudotime values are converted into "effectorness" measurements, which are comparable across cytokine conditions

# LOADING R LIBRARIES
library(Seurat)
library(monocle)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(VennDiagram)
library(reshape2)

# LOADING DATA
## Loading a Seurat object containing the expression matrix and annotations for all cells in the study
## Please refer to "scRNAseq_exploratory_data_analysis.R" to see how this object was generated
T.cells <- readRDS("allCells_Seurat.rds")

## Loading a list of cell cycle genes define by the Regev lab
cc.genes <- readLines("~/cytoimmgen/phase_1/data/regev_lab_cell_cycle_genes.txt")

# SPLITTING DATA BY CYTOKINE CONDITION
uns.cells <- SubsetData(T.cells, cells.use=T.cells@cell.names[T.cells@meta.data$cytokine.condition == "UNS"])
Th0 <- SubsetData(T.cells, cells.use = T.cells@cell.names[T.cells@meta.data$cytokine.condition=="Th0"])
Th2 <- SubsetData(T.cells, cells.use = T.cells@cell.names[T.cells@meta.data$cytokine.condition=="Th2"])
Th17 <- SubsetData(T.cells, cells.use = T.cells@cell.names[T.cells@meta.data$cytokine.condition=="Th17"])
iTreg <- SubsetData(T.cells, cells.use = T.cells@cell.names[T.cells@meta.data$cytokine.condition=="iTreg"])

# IDENTIFYING HIGHLY VARIABLE GENES
uns.cells <- ScaleData(uns.cells, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))
Th0 <- ScaleData(Th0, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))
Th2 <- ScaleData(Th2, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))
Th17 <- ScaleData(Th17, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))
iTreg <- ScaleData(iTreg, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))

uns.cells <- FindVariableGenes(uns.cells,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)
Th0 <- FindVariableGenes(Th0,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)
Th2 <- FindVariableGenes(Th2,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)
Th17 <- FindVariableGenes(Th17,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)
iTreg <- FindVariableGenes(iTreg,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)

# BUILDING A PSEUDOTIME TRAJECTORY FOR EACH CYTOKINE CONDITION
## Converting objects to the CellDataSet format required by Monocle 2
UNS.monocle <- importCDS(uns.cells)
Th0.monocle <- importCDS(Th0)
Th2.monocle <- importCDS(Th2)
Th17.monocle <- importCDS(Th17)
iTreg.monocle <- importCDS(iTreg)

## Estimating size factors
UNS.monocle <- estimateSizeFactors(UNS.monocle)
UNS.monocle <- estimateDispersions(UNS.monocle)

Th0.monocle <- estimateSizeFactors(Th0.monocle)
Th0.monocle <- estimateDispersions(Th0.monocle)

Th2.monocle <- estimateSizeFactors(Th2.monocle)
Th2.monocle <- estimateDispersions(Th2.monocle)

Th17.monocle <- estimateSizeFactors(Th17.monocle)
Th17.monocle <- estimateDispersions(Th17.monocle)

iTreg.monocle <- estimateSizeFactors(iTreg.monocle)
iTreg.monocle <- estimateDispersions(iTreg.monocle)

## Using highly variable genes defined by Seurat to order cells in trajectories.
### Any genes involved in cell cycle are removed.
ordering_genes.UNS <- uns.cells@var.genes
ordering_genes.UNS <- setdiff(ordering_genes.UNS, unique(unlist(cc.genes)))
UNS.monocle <- setOrderingFilter(UNS.monocle, ordering_genes.UNS)

ordering_genes <- Th0@var.genes
ordering_genes <- setdiff(ordering_genes, unique(unlist(cc.genes)))
Th0.monocle <- setOrderingFilter(Th0.monocle, ordering_genes)

ordering_genes.Th2 <- Th2@var.genes
ordering_genes.Th2 <- setdiff(ordering_genes.Th2, unique(unlist(cc.genes)))
Th2.monocle <- setOrderingFilter(Th2.monocle, ordering_genes.Th2)

ordering_genes.Th17 <- Th17@var.genes
ordering_genes.Th17 <- setdiff(ordering_genes.Th17, unique(unlist(cc.genes)))
Th17.monocle <- setOrderingFilter(Th17.monocle, ordering_genes.Th17)

ordering_genes.iTreg <- iTreg@var.genes
ordering_genes.iTreg <- setdiff(ordering_genes.iTreg, unique(unlist(cc.genes)))
iTreg.monocle <- setOrderingFilter(iTreg.monocle, ordering_genes.iTreg)

## Reducing dimensions with the DDRTree algorithm
UNS.monocle <- reduceDimension(UNS.monocle, max_components = 2,method = 'DDRTree')
Th0.monocle <- reduceDimension(Th0.monocle, max_components = 2,method = 'DDRTree')
Th2.monocle <- reduceDimension(Th2.monocle, max_components = 2,method = 'DDRTree')
Th17.monocle <- reduceDimension(Th17.monocle, max_components = 2,method = 'DDRTree')
iTreg.monocle <- reduceDimension(iTreg.monocle, max_components = 2,method = 'DDRTree')

## Ordering cells into pseudotime trajectories
UNS.monocle <- orderCells(UNS.monocle)
Th0.monocle <- orderCells(Th0.monocle)
Th2.monocle <- orderCells(Th2.monocle)
Th17.monocle <- orderCells(Th17.monocle)
iTreg.monocle <- orderCells(iTreg.monocle)

### Reversing trajectory orders whenever necessary (as determined by visualisation)
Th2.monocle <- orderCells(Th2.monocle, reverse=T)
Th17.monocle <- orderCells(Th17.monocle, reverse=T)
iTreg.monocle <- orderCells(iTreg.monocle, reverse=T)

## Identifying genes associated with pseudotime in each trajectory
diff_test_res.UNS <- differentialGeneTest(UNS.monocle[ordering_genes.UNS,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names.UNS <- row.names(subset(diff_test_res.UNS, qval < 0.01))

diff_test_res <- differentialGeneTest(Th0.monocle[ordering_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))

diff_test_res.Th2 <- differentialGeneTest(Th2.monocle[ordering_genes.Th2,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names.Th2 <- row.names(subset(diff_test_res.Th2, qval < 0.01))

diff_test_res.Th17 <- differentialGeneTest(Th17.monocle[ordering_genes.Th17,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names.Th17 <- row.names(subset(diff_test_res.Th17, qval < 0.01))

diff_test_res.iTreg <- differentialGeneTest(iTreg.monocle[ordering_genes.iTreg,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names.iTreg <- row.names(subset(diff_test_res.iTreg, qval < 0.01))

# CONVERTING PSEUDOTIME INTO EFFECTORNESS
## Definig effectorness as the pseudotime values of each cel scaled to the [0,1] range. 
## This make values comparable across cytokine conditions
UNS.Eff <- (pData(UNS.monocle)$Pseudotime - min(pData(UNS.monocle)$Pseudotime))/(max(pData(UNS.monocle)$Pseudotime)-min(pData(UNS.monocle)$Pseudotime))
Th0.Eff <- (pData(Th0.monocle)$Pseudotime - min(pData(Th0.monocle)$Pseudotime))/(max(pData(Th0.monocle)$Pseudotime)-min(pData(Th0.monocle)$Pseudotime))
Th2.Eff <- (pData(Th2.monocle)$Pseudotime - min(pData(Th2.monocle)$Pseudotime))/(max(pData(Th2.monocle)$Pseudotime)-min(pData(Th2.monocle)$Pseudotime))
Th17.Eff <- (pData(Th17.monocle)$Pseudotime - min(pData(Th17.monocle)$Pseudotime))/(max(pData(Th17.monocle)$Pseudotime)-min(pData(Th17.monocle)$Pseudotime))
iTreg.Eff <- (pData(iTreg.monocle)$Pseudotime - min(pData(iTreg.monocle)$Pseudotime))/(max(pData(iTreg.monocle)$Pseudotime)-min(pData(iTreg.monocle)$Pseudotime))

names(UNS.Eff) <- uns.cells@cell.names
names(Th0.Eff) <- Th0@cell.names
names(Th2.Eff) <- Th2@cell.names
names(Th17.Eff) <- Th17@cell.names
names(iTreg.Eff) <- iTreg@cell.names

## Adding the estimated effectorness values as metadata to the original Seurat object
### Adding an "effectorness" variable to the metadta and setting it to a dummy value
T.cells@meta.data$effectorness <- NA

### Adding the correct "effectorness" value to each cell in each cytokine condition
T.cells@meta.data[names(UNS.Eff),]$effectorness <- UNS.Eff
T.cells@meta.data[names(Th0.Eff),]$effectorness <- Th0.Eff
T.cells@meta.data[names(Th2.Eff),]$effectorness <- Th2.Eff
T.cells@meta.data[names(Th17.Eff),]$effectorness <- Th17.Eff
T.cells@meta.data[names(iTreg.Eff),]$effectorness <- iTreg.Eff

# EXPORTING RESULTS
saveRDS(T.cells, "allCells_Seurat_with_effectorness.rds")
