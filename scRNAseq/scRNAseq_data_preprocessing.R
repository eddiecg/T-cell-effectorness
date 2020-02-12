### PRE-PROCESSING SINGLE-CELL RNA-SEQ DATA AND GENERATING A CELL ANNOTATIONS METADTA TABLE ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code shows how 10X count tables (after cellranger analysis) were used to generate the matrix.mtx, genes.tsv, barcodes.tsv, and metadata.txt files
# These four files are available for download as supplementary material of the study "Single-cell transcriptomics identifies an effectorness gradient shaping the response of CD4+ T cells to cytokines"

# LOADING R LIBRARIES
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)

# LOADING 10X DATA
## Loading expression matrices for each 10X sample in the study
##Each sample corresponds to a replicate of one cytokine condition
### Naive T cells
N.UNS.data <- Read10X("./cellranger/I1465/filtered_gene_bc_matrices/GRCh38/")
N.Th0.data <- Read10X("./cellranger/I1456/filtered_gene_bc_matrices/GRCh38/")
N.Th2.data <- Read10X("./cellranger/I1458/filtered_gene_bc_matrices/GRCh38/")
N.Th17.data.1 <- Read10X("./cellranger/I1462/filtered_gene_bc_matrices/GRCh38/")
N.Th17.data.2 <- Read10X("./cellranger/I1470/filtered_gene_bc_matrices/GRCh38/")
N.iTreg.data.1 <- Read10X("./cellranger/I1460/filtered_gene_bc_matrices/GRCh38/")
N.iTreg.data.2 <- Read10X("./cellranger/I1468/filtered_gene_bc_matrices/GRCh38/")

### Memory T cells
M.UNS.data.1 <- Read10X("./cellranger/I1466/filtered_gene_bc_matrices/GRCh38/")
M.UNS.data.2 <- Read10X("./cellranger/I1467/filtered_gene_bc_matrices/GRCh38/")
M.Th0.data.1 <- Read10X("./cellranger/I1457/filtered_gene_bc_matrices/GRCh38/")
M.Th0.data.2 <- Read10X("./cellranger/I1464/filtered_gene_bc_matrices/GRCh38/")
M.Th2.data <- Read10X("./cellranger/I1459/filtered_gene_bc_matrices/GRCh38/")
M.Th17.data.1 <- Read10X("./cellranger/I1463/filtered_gene_bc_matrices/GRCh38/")
M.Th17.data.2 <- Read10X("./cellranger/I1471/filtered_gene_bc_matrices/GRCh38/")
M.iTreg.data.1 <- Read10X("./cellranger/I1461/filtered_gene_bc_matrices/GRCh38/")
M.iTreg.data.2 <- Read10X("./cellranger/I1469/filtered_gene_bc_matrices/GRCh38/")

# CREATING SEURAT OBJECTS
## Naive T cells
N.UNS <- CreateSeuratObject(N.UNS.data, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
N.Th0 <- CreateSeuratObject(N.Th0.data, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
N.Th2 <- CreateSeuratObject(N.Th2.data, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
N.Th17.1 <- CreateSeuratObject(N.Th17.data.1, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
N.Th17.2 <- CreateSeuratObject(N.Th17.data.2, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
N.iTreg.1 <- CreateSeuratObject(N.iTreg.data.1, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
N.iTreg.2 <- CreateSeuratObject(N.iTreg.data.2, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")

## Memory T cells
M.UNS.1 <- CreateSeuratObject(M.UNS.data.1, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.UNS.2 <- CreateSeuratObject(M.UNS.data.2, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.Th0.1 <- CreateSeuratObject(M.Th0.data.1, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.Th0.2 <- CreateSeuratObject(M.Th0.data.2, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.Th2 <- CreateSeuratObject(M.Th2.data, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.Th17.1 <- CreateSeuratObject(M.Th17.data.1, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.Th17.2 <- CreateSeuratObject(M.Th17.data.2, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.iTreg.1 <- CreateSeuratObject(M.iTreg.data.1, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")
M.iTreg.2 <- CreateSeuratObject(M.iTreg.data.2, min.cells = 3, min.genes = 200, project = "T_cell_effectorness")

# ASSIGNING CELLS TO INDIVIDUALS
## Loading genotype deconvolution results from Cardelino
## Please refer to "scRNAseq_cardelino_get_donor_ids.R" to see how these files were generated
### Naive T cells
N.UNS.donors <- readRDS("./donor_deconvolution/I1465.donorID.rds")
N.Th0.donors <- readRDS("./donor_deconvolution/I1456.donorID.rds")
N.Th2.donors <- readRDS("./donor_deconvolution/I1458.donorID.rds")
N.Th17.donors.1 <- readRDS("./donor_deconvolution/I1462.donorID.rds")
N.Th17.donors.2 <- readRDS("./donor_deconvolution/I1470.donorID.rds")
N.iTreg.donors.1 <- readRDS("./donor_deconvolution/I1460.donorID.rds")
N.iTreg.donors.2 <- readRDS("./donor_deconvolution/I1468.donorID.rds")

### Memory T cells
M.UNS.donors.1 <- readRDS("./donor_deconvolution/I1466.donorID.rds")
M.UNS.donors.2 <- readRDS("./donor_deconvolution/I1467.donorID.rds")
M.Th0.donors.1 <- readRDS("./donor_deconvolution/I1457.donorID.rds")
M.Th0.donors.2 <- readRDS("./donor_deconvolution/I1464.donorID.rds")
M.Th2.donors <- readRDS("./donor_deconvolution/I1459.donorID.rds")
M.Th17.donors.1 <- readRDS("./donor_deconvolution/I1463.donorID.rds")
M.Th17.donors.2 <- readRDS("./donor_deconvolution/I1471.donorID.rds")
M.iTreg.donors.1 <- readRDS("./donor_deconvolution/I1461.donorID.rds")
M.iTreg.donors.2 <- readRDS("./donor_deconvolution/I1469.donorID.rds")

## Loading a list of all donor IDs in the study
donor.ids <- read.csv("donor_ids.csv")

## Defining a function which matches the donor assigned to each cell to its donor ID in the study
getDonorIDs <- function(assigned_cells, donor_ids, sample_id){
  d <- assigned_cells$donor_id
  sapply(assigned_cells$donor_id, FUN=function(d){
    donor=ifelse(d=="unassigned", "unassigned",
                 ifelse(d=="doublet", "doublet",
                        as.character(donor_ids[donor_ids$sample_id==sample_id & donor_ids$donor==d,]$donor_id)))
    return(donor)
  })
}

## Assigning donor IDs to cells
### Naive T cells
N.UNS@meta.data$donor.id <- getDonorIDs(N.UNS.donors$assigned, donor.ids, "I1465")
N.Th0@meta.data$donor.id <- getDonorIDs(N.Th0.donors$assigned, donor.ids, "I1456")
N.Th2@meta.data$donor.id <- getDonorIDs(N.Th2.donors$assigned, donor.ids, "I1458")
N.Th17.1@meta.data$donor.id <- getDonorIDs(N.Th17.donors.1$assigned, donor.ids, "I1462")
N.Th17.2@meta.data$donor.id <- getDonorIDs(N.Th17.donors.2$assigned, donor.ids, "I1470")
N.iTreg.1@meta.data$donor.id <- getDonorIDs(N.iTreg.donors.1$assigned, donor.ids, "I1460")
N.iTreg.2@meta.data$donor.id <- getDonorIDs(N.iTreg.donors.2$assigned, donor.ids, "I1468")

### Memory T cells
M.UNS.1@meta.data$donor.id <- getDonorIDs(M.UNS.donors.1$assigned, donor.ids, "I1467")
M.UNS.2@meta.data$donor.id <- getDonorIDs(M.UNS.donors.2$assigned, donor.ids, "I1466")
M.Th0.1@meta.data$donor.id <- getDonorIDs(M.Th0.donors.1$assigned, donor.ids, "I1457")
M.Th0.2@meta.data$donor.id <- getDonorIDs(M.Th0.donors.2$assigned, donor.ids, "I1464")
M.Th2@meta.data$donor.id <- getDonorIDs(M.Th2.donors$assigned, donor.ids, "I1459")
M.Th17.1@meta.data$donor.id <- getDonorIDs(M.Th17.donors.1$assigned, donor.ids, "I1463")
M.Th17.2@meta.data$donor.id <- getDonorIDs(M.Th17.donors.2$assigned, donor.ids, "I1471")
M.iTreg.1@meta.data$donor.id <- getDonorIDs(M.iTreg.donors.1$assigned, donor.ids, "I1461")
M.iTreg.2@meta.data$donor.id <- getDonorIDs(M.iTreg.donors.2$assigned, donor.ids, "I1469")

# MERGING TECHNICAL REPLICATES OF THE SAME CONDITIONS
N.Th17 <- MergeSeurat(N.Th17.1, N.Th17.2, add.cell.id1 = "r1", add.cell.id2 = "r2", project = "Cytoimmgen_10X_optimization")
N.iTreg <- MergeSeurat(N.iTreg.1, N.iTreg.2, add.cell.id1 = "r1", add.cell.id2 = "r2", project = "Cytoimmgen_10X_optimization")
M.UNS <- MergeSeurat(M.UNS.1, M.UNS.2, add.cell.id1 = "r1", add.cell.id2 = "r2", project = "Cytoimmgen_10X_optimization")
M.Th0 <- MergeSeurat(M.Th0.1, M.Th0.2, add.cell.id1 = "r1", add.cell.id2 = "r2", project = "Cytoimmgen_10X_optimization")
M.Th17 <- MergeSeurat(M.Th17.1, M.Th17.2, add.cell.id1 = "r1", add.cell.id2 = "r2", project = "Cytoimmgen_10X_optimization")
M.iTreg <- MergeSeurat(M.iTreg.1, M.iTreg.2, add.cell.id1 = "r1", add.cell.id2 = "r2", project = "Cytoimmgen_10X_optimization")

### Removing cells from donor 'D2' from all the iTreg samples in the study (due to unusually high dropout rate of this donot in this condition)
keep <- rownames(N.iTreg@meta.data[N.iTreg@meta.data$donor %in% c("D1","D3","D4"),])
N.iTreg <- SubsetData(N.iTreg, cells.use = keep)
rm(keep)

# CREATING METADATA FIELDS
## Batch number
N.UNS@meta.data$batch_10X <- c(rep(2,ncol(N.UNS@raw.data)))
N.Th0@meta.data$batch_10X <- c(rep(1,ncol(N.Th0@raw.data)))
N.Th2@meta.data$batch_10X <- c(rep(1,ncol(N.Th2@raw.data)))
N.Th17@meta.data$batch_10X <- c(rep(1,ncol(N.Th17.1@raw.data)),rep(2,ncol(N.Th17.2@raw.data)))
N.iTreg@meta.data$batch_10X <- c(rep(1,ncol(N.iTreg.1@raw.data)),rep(2,ncol(N.iTreg.2@raw.data)))

M.UNS@meta.data$batch_10X <- c(rep(2,ncol(M.UNS.1@raw.data)),rep(2,ncol(M.UNS.2@raw.data)))
M.Th0@meta.data$batch_10X <- c(rep(1,ncol(M.Th0.1@raw.data)),rep(2,ncol(M.Th0.2@raw.data)))
M.Th2@meta.data$batch_10X <- c(rep(1,ncol(M.Th2@raw.data)))
M.Th17@meta.data$batch_10X <- c(rep(1,ncol(M.Th17.1@raw.data)),rep(2,ncol(M.Th17.2@raw.data)))
M.iTreg@meta.data$batch_10X <- c(rep(1,ncol(M.iTreg.1@raw.data)),rep(2,ncol(M.iTreg.2@raw.data)))

## Cell type
N.UNS@meta.data$cell.type <- c(rep("Naive",ncol(N.UNS@raw.data)))
N.Th0@meta.data$cell.type <- c(rep("Naive",ncol(N.Th0@raw.data)))
N.Th2@meta.data$cell.type <- c(rep("Naive",ncol(N.Th2@raw.data)))
N.Th17@meta.data$cell.type <- c(rep("Naive",ncol(N.Th17@raw.data)))
N.iTreg@meta.data$cell.type <- c(rep("Naive",ncol(N.iTreg@raw.data)))

M.UNS@meta.data$cell.type <- c(rep("Memory",ncol(M.UNS@raw.data)))
M.Th0@meta.data$cell.type <- c(rep("Memory",ncol(M.Th0@raw.data)))
M.Th2@meta.data$cell.type <- c(rep("Memory",ncol(M.Th2@raw.data)))
M.Th17@meta.data$cell.type <- c(rep("Memory",ncol(M.Th17@raw.data)))
M.iTreg@meta.data$cell.type <- c(rep("Memory",ncol(M.iTreg@raw.data)))

## Cytokine condition
N.UNS@meta.data$cytokine.condition <- c(rep("UNS",ncol(N.UNS@raw.data)))
N.Th0@meta.data$cytokine.condition <- c(rep("Th0",ncol(N.Th0@raw.data)))
N.Th2@meta.data$cytokine.condition <- c(rep("Th2",ncol(N.Th2@raw.data)))
N.Th17@meta.data$cytokine.condition <- c(rep("Th17",ncol(N.Th17@raw.data)))
N.iTreg@meta.data$cytokine.condition <- c(rep("iTreg",ncol(N.iTreg@raw.data)))

M.UNS@meta.data$cytokine.condition <- c(rep("UNS",ncol(M.UNS@raw.data)))
M.Th0@meta.data$cytokine.condition <- c(rep("Th0",ncol(M.Th0@raw.data)))
M.Th2@meta.data$cytokine.condition <- c(rep("Th2",ncol(M.Th2@raw.data)))
M.Th17@meta.data$cytokine.condition <- c(rep("Th17",ncol(M.Th17@raw.data)))
M.iTreg@meta.data$cytokine.condition <- c(rep("iTreg",ncol(M.iTreg@raw.data)))

## Mitochondrial content
addMito <- function(dat){
  mito.genes <- grep(pattern = "^MT-", x=rownames(dat@data), value=TRUE)
  percent.mito <- colSums(as.matrix(dat@raw.data[mito.genes,]))/colSums(as.matrix(dat@raw.data))
  dat <- AddMetaData(dat, metadata = percent.mito, col.name = "percent.mito")
  return(dat)
}

N.UNS <- addMito(N.UNS)
N.Th0 <- addMito(N.Th0)
N.Th2 <- addMito(N.Th2)
N.Th17 <- addMito(N.Th17)
N.iTreg <- addMito(N.iTreg)

M.UNS <- addMito(M.UNS)
M.Th0 <- addMito(M.Th0)
M.Th2 <- addMito(M.Th2)
M.Th17 <- addMito(M.Th17)
M.iTreg <- addMito(M.iTreg)

# FILTERING OUT LOW QUALITY CELLS CELLS
## Removing doublets and unassigned cells (based on the results from Cardelino)
removeUnassigned <- function(dat){
  assigned <- rownames(dat@meta.data[dat@meta.data$donor.id %in% c("D1","D2","D3","D4"),])
  dat <- SubsetData(dat, cells.use = assigned)
  return(dat)
}

N.UNS <- removeUnassigned(N.UNS)
N.Th0 <- removeUnassigned(N.Th0)
N.Th2 <- removeUnassigned(N.Th2)
N.Th17 <- removeUnassigned(N.Th17)
N.iTreg <- removeUnassigned(N.iTreg)

M.UNS <- removeUnassigned(M.UNS)
M.Th0 <- removeUnassigned(M.Th0)
M.Th2 <- removeUnassigned(M.Th2)
M.Th17 <- removeUnassigned(M.Th17)
M.iTreg <- removeUnassigned(M.iTreg)

## Filtering cells by mitochondrial (<7.5%) and UMI content (>500)
N.UNS <- FilterCells(N.UNS, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
N.Th0 <- FilterCells(N.Th0, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
N.Th2 <- FilterCells(N.Th2, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
N.Th17 <- FilterCells(N.Th17, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
N.iTreg <- FilterCells(N.iTreg, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))

M.UNS <- FilterCells(M.UNS, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
M.Th0 <- FilterCells(M.Th0, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
M.Th2 <- FilterCells(M.Th2, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
M.Th17 <- FilterCells(M.Th17, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))
M.iTreg <- FilterCells(M.iTreg, subset.names = c("nGene","percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.075))

# MERGING DATA SETS FROM ALL CONDITIONS INTO A SINGLE OBJECT
cells <- MergeSeurat(N.UNS, N.Th0, add.cell.id1=paste("N",N.UNS@meta.data$cytokine.condition,sep="_"), add.cell.id2=paste("N",N.Th0@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, N.Th2, add.cell.id1=NULL, add.cell.id2=paste("N",N.Th2@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, N.Th17, add.cell.id1=NULL, add.cell.id2=paste("N",N.Th17@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, N.iTreg, add.cell.id1=NULL, add.cell.id2=paste("N",N.iTreg@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, M.UNS, add.cell.id1=NULL, add.cell.id2=paste("M",M.UNS@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, M.Th0, add.cell.id1=NULL, add.cell.id2=paste("M",M.Th0@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, M.Th2, add.cell.id1=NULL, add.cell.id2=paste("M",M.Th2@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, M.Th17, add.cell.id1=NULL, add.cell.id2=paste("M",M.Th17@meta.data$cytokine.condition,sep="_"))
cells <- MergeSeurat(cells, M.iTreg, add.cell.id1=NULL, add.cell.id2=paste("M",M.iTreg@meta.data$cytokine.condition,sep="_"))

# ASSIGNING A CELL CYCLE SCORE TO EACH CELL
## Loading a list of cell cycle genes compiled by the Regev lab 
cc.genes <- readLines("regev_lab_cell_cycle_genes.txt")
S.genes <- cc.genes[1:43]
G2M.genes <- cc.genes[44:97]

## Scoring genes by cell cycle gene expression
all.cells <- NormalizeData(all.cells, normalization.method = "LogNormalize")
all.cells <- CellCycleScoring(all.cells, s.genes=S.genes, g2m.genes=G2M.genes, set.ident=TRUE)

# EXPORTING FINAL, HIGH QUALITY DATA SET TO .MTX AND .TSV FILES
writeMM(cells@raw.data, "NCOMMS-19-7936188_scRNAseq_raw_UMIs.mtx")
write.table(rownames(cells@raw.data), "NCOMMS-19-7936188_scRNAseq_genes.tsv",quote=F,sep="\t", row.names=F, col.names=F)
write.table(cells@cell.names, "NCOMMS-19-7936188_scRNAseq_barcodes.tsv",quote=F,sep="\t", row.names=F, col.names=F)
write.table(cells@meta.data, "NCOMMS-19-7936188_scRNAseq_metadata.txt",quote=F,sep="\t")

### NOTE: The final files available for download contain additional metadta fields (from further analysis steps such as clustering or trajectory analysis).
###       Futhermore, their row names were slightly modified to make the data more readable
###       A small proportion of contaminating erythroid cells was also removed from the final version of these files

