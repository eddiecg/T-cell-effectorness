### INFERRING CELL STATE-SPECIFIC SIGNATURES FROM MATCHING RNA AND PROTEIN DATA ###
## Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)

# INSTALLING PROTEOGENOMIC LIBRARY
## Run the following line to install the "proteogenomic" R package
library(devtools)
install_github("eddiecg/proteogenomic")

# LOADING LIBRARIES
library(DESeq2)
library(proteogenomic)

# LOADING DATA
## Loading lists of differentially expressed gene-protein pairs (fCI)
### Naive T cells
Th1.N <- read.csv(file="naive_Th1_vs_Th0_fCI_results.csv")
Th2.N <- read.csv(file="naive_Th2_vs_Th0_fCI_results.csv")
Th17.N <- read.csv(file="naive_Th17_vs_Th0_fCI_results.csv")
iTreg.N <- read.csv(file="naive_iTreg_vs_Th0_fCI_results.csv")
IFNb.N <- read.csv(file="naive_IFNB_vs_Th0_fCI_results.csv")

### Memory T cells
Th1.M <- read.csv(file="memory_Th1_vs_Th0_fCI_results.csv")
Th2.M <- read.csv(file="memory_Th2_vs_Th0_fCI_results.csv")
Th17.M <- read.csv(file="memory_Th17_vs_Th0_fCI_results.csv")
iTreg.M <- read.csv(file="memory_iTreg_vs_Th0_fCI_results.csv")
IFNb.M <- read.csv(file="memory_IFNB_vs_Th0_fCI_results.csv")

## Merging all condition-specific differentially expressed genes into a single list per cell type
DE_genes_naive <- unique(c(as.character(Th1.N$DEG_Names),
                           as.character(Th2.N$DEG_Names),
                           as.character(Th17.N$DEG_Names),
                           as.character(iTreg.N$DEG_Names),
                           as.character(IFNb.N$DEG_Names))
                         )

DE_genes_memory <- unique(c(as.character(Th1.M$DEG_Names),
                            as.character(Th2.M$DEG_Names),
                            as.character(Th17.M$DEG_Names),
                            as.character(iTreg.M$DEG_Names),
                            as.character(IFNb.M$DEG_Names))
                          )
					
## Loading RNA-seq and proteomics count tables
### Loading RNA counts
dds <- readRDS("rawCounts_bulkRNAseq_DESeq2.rds")

### Loading protein abundances
prots <- readRDS("proteinAbundances_summarizedExperiment.rds")

## Pre-processing RNA-seq and proteomics data
### Subsetting data by cell type and time point 
dds.naive <- dds[,(dds$cell_type == "CD4_Naive") & (dds$stimulation_time == "5d")]
dds.memory <- dds[,(dds$cell_type == "CD4_Memory") & (dds$stimulation_time == "5d")]

prots.N <- prots[,prots$cell_type == "CD4_naive" ]
prots.M <- prots[,prots$cell_type == "CD4_memory" ]

#### Subsetting RNA-seq data to differentially expressed genes (fCI) only
rna.N <- dds.naive[subset(rowData(dds.naive), Gene_symbol %in% DE_genes_naive)$Gene_id,]
rownames(rna.N) <- rowData(rna.N)$Gene_symbol

rna.M <- dds.memory[subset(rowData(dds.memory), Gene_symbol %in% DE_genes_memory)$Gene_id,]
rownames(rna.M) <- rowData(rna.M)$Gene_symbol

#### Subsetting proteomics data to differentially expressed genes (fCI) only
prots.N <- prots.N[as.character(subset(rowData(prots.N), Gene_name %in% DE_genes_naive)$Gene_name),]
rownames(prots.N) <- rowData(prots.N)$Gene_name

prots.M <- prots.M[as.character(subset(rowData(prots.M), Gene_name %in% DE_genes_memory)$Gene_name),]
rownames(prots.M) <- rowData(prots.M)$Gene_name

#### Matching RNA and protein data sets
prots.N <- prots.N[intersect(rownames(rna.N), rownames(prots.N)),]
prots.M <- prots.M[intersect(rownames(rna.M),rownames(prots.M)),]

#### Removing duplicates
rna.N <- rna.N[!duplicated(rownames(rna.N)),]
rna.M <- rna.M[!duplicated(rownames(rna.M)),]

# FORMATTING DATA SETS AND MATCHING COLUMNS
### Extracting counts/abundances
prots.N <- data.frame(assay(prots.N))
rna.N <- counts(rna.N, normalized=T)

prots.M <- data.frame(assay(prots.M))
rna.M <- counts(rna.M, normalized=T)

### Matching column order between data sets
rna.N <- rna.N[,c(grep("Th0",colnames(rna.N)),
                  grep("Th1",colnames(rna.N)),
                  grep("Th2",colnames(rna.N)),
                  grep("Th17",colnames(rna.N)),
                  grep("iTreg",colnames(rna.N)),
                  grep("IFNB",colnames(rna.N))
                  )]
prots.N <- prots.N[,c(grep("Th0",colnames(prots.N)),
                      grep("Th1",colnames(prots.N)),
                      grep("Th2",colnames(prots.N)),
                      grep("Th17",colnames(prots.N)),
                      grep("iTreg",colnames(prots.N)),
                      grep("IFNB",colnames(prots.N))
                      )]

rna.M <- rna.M[,c(grep("Th0",colnames(rna.M)),
                  grep("Th1",colnames(rna.M)),
                  grep("Th2",colnames(rna.M)),
                  grep("Th17",colnames(rna.M)),
                  grep("iTreg",colnames(rna.M)),
                  grep("IFNB",colnames(rna.M))
                  )]
prots.M <- prots.M[,c(grep("Th0",colnames(prots.M)),
                      grep("Th1",colnames(prots.M)),
                      grep("Th2",colnames(prots.M)),
                      grep("Th17",colnames(prots.M)),
                      grep("iTreg",colnames(prots.M)),
                      grep("IFNB",colnames(prots.M))
                      )]

### Trimming column names
colnames(prots.N) <- gsub("\\.1","",colnames(prots.N))
colnames(prots.M) <- gsub("\\.1","",colnames(prots.M))

### Simplifying column names to only be the cytokine condition
colnames(rna.N) <-  as.character(sapply(colnames(rna.N), FUN=function(x){strsplit(x, split="_")[[1]][5]}))
colnames(prots.N) <- as.character(sapply(colnames(prots.N), FUN=function(x){strsplit(x, split="_")[[1]][4]}))

colnames(rna.M) <-  as.character(sapply(colnames(rna.M), FUN=function(x){strsplit(x, split="_")[[1]][5]}))
colnames(prots.M) <- as.character(sapply(colnames(prots.M), FUN=function(x){strsplit(x, split="_")[[1]][4]}))

### Adding "dummy columns" (i.e. columns with NA values for every gene) for all missing samples
#### Naive IFNB (RNA)
rna.N <- data.frame(rna.N,IFNB=rep(NA,nrow(rna.N)))
colnames(rna.N) <- gsub("\\..*","",colnames(rna.N))

#### Memory IFNB (proteomics)
prots.M <- data.frame(prots.M,IFNB=rep(NA,nrow(prots.M)))
colnames(prots.M) <- gsub("\\..*","",colnames(prots.M))


## MIXIMG TH17- AND ITREG-STIMULATED MEMORY CELLS INTO A SINGLE CATEGORY
### RNA
colnames(rna.M) <- gsub("Mixed","Th17_iTreg",gsub("iTreg","Mixed",gsub("Th17","Mixed",colnames(rna.M))))
colnames(prots.M) <- gsub("Mixed","Th17_iTreg",gsub("iTreg","Mixed",gsub("Th17","Mixed",colnames(prots.M))))

# CALCULATING SPECIFICITIES
## Naive T cells
N.specs <- testSpecificities(rna.N, prots.N, sample.labels=c("Th1","Th2","Th17","iTreg","IFNB"),iter=10000)

## Memory T cells
M.specs <- testSpecificities(rna.M, prots.M, sample.labels=c("Th1","Th17_iTreg","IFNB"),iter=10000)

# RETRIEVEING CELL STATE SPECIFIC SIGNATURES
fdr.thres <- 0.1
sp.thres <- 0.7

## Naive T cells
Th1.N.signature <- data.frame(gene=rownames(N.specs$specificities[N.specs$p.adj[,"Th1"] < fdr.Thres & N.specs$specificities[,"Th1"] > sp.Thres,]),
                              specificity=N.specs$specificities[N.specs$p.adj[,"Th1"] < fdr.Thres & N.specs$specificities[,"Th1"] > sp.Thres,]$Th1,
                              p.value=N.specs$p.val[N.specs$p.adj[,"Th1"] < fdr.Thres & N.specs$specificities[,"Th1"] > sp.Thres,]$Th1,
                              p.adj=N.specs$p.adj[N.specs$p.adj[,"Th1"] < fdr.Thres & N.specs$specificities[,"Th1"] > sp.Thres,]$Th1)

Th2.N.signature <- data.frame(gene=rownames(N.specs$specificities[N.specs$p.adj[,"Th2"] < fdr.Thres & N.specs$specificities[,"Th2"] > sp.Thres,]),
                              specificity=N.specs$specificities[N.specs$p.adj[,"Th2"] < fdr.Thres & N.specs$specificities[,"Th2"] > sp.Thres,]$Th2,
                              p.value=N.specs$p.val[N.specs$p.adj[,"Th2"] < fdr.Thres & N.specs$specificities[,"Th2"] > sp.Thres,]$Th2,
                              p.adj=N.specs$p.adj[N.specs$p.adj[,"Th2"] < fdr.Thres & N.specs$specificities[,"Th2"] > sp.Thres,]$Th2)

Th17.N.signature <- data.frame(gene=rownames(N.specs$specificities[N.specs$p.adj[,"Th17"] < fdr.Thres & N.specs$specificities[,"Th17"] > sp.Thres,]),
                              specificity=N.specs$specificities[N.specs$p.adj[,"Th17"] < fdr.Thres & N.specs$specificities[,"Th17"] > sp.Thres,]$Th17,
                              p.value=N.specs$p.val[N.specs$p.adj[,"Th17"] < fdr.Thres & N.specs$specificities[,"Th17"] > sp.Thres,]$Th17,
                              p.adj=N.specs$p.adj[N.specs$p.adj[,"Th17"] < fdr.Thres & N.specs$specificities[,"Th17"] > sp.Thres,]$Th17)

iTreg.N.signature <- data.frame(gene=rownames(N.specs$specificities[N.specs$p.adj[,"iTreg"] < fdr.Thres & N.specs$specificities[,"iTreg"] > sp.Thres,]),
                              specificity=N.specs$specificities[N.specs$p.adj[,"iTreg"] < fdr.Thres & N.specs$specificities[,"iTreg"] > sp.Thres,]$iTreg,
                              p.value=N.specs$p.val[N.specs$p.adj[,"iTreg"] < fdr.Thres & N.specs$specificities[,"iTreg"] > sp.Thres,]$iTreg,
                              p.adj=N.specs$p.adj[N.specs$p.adj[,"iTreg"] < fdr.Thres & N.specs$specificities[,"iTreg"] > sp.Thres,]$iTreg)

IFNB.N.signature <- data.frame(gene=rownames(N.specs$specificities[N.specs$p.adj[,"IFNB"] < fdr.Thres & N.specs$specificities[,"IFNB"] > sp.Thres,]),
                                specificity=N.specs$specificities[N.specs$p.adj[,"IFNB"] < fdr.Thres & N.specs$specificities[,"IFNB"] > sp.Thres,]$IFNB,
                                p.value=N.specs$p.val[N.specs$p.adj[,"IFNB"] < fdr.Thres & N.specs$specificities[,"IFNB"] > sp.Thres,]$IFNB,
                                p.adj=N.specs$p.adj[N.specs$p.adj[,"IFNB"] < fdr.Thres & N.specs$specificities[,"IFNB"] > sp.Thres,]$IFNB)


Th1.N.signature <- Th1.N.signature[order(-Th1.N.signature$specificity),]
Th2.N.signature <- Th2.N.signature[order(-Th2.N.signature$specificity),]
Th17.N.signature <- Th17.N.signature[order(-Th17.N.signature$specificity),]
iTreg.N.signature <- iTreg.N.signature[order(-iTreg.N.signature$specificity),]
IFNB.N.signature <- IFNB.N.signature[order(-IFNB.N.signature$specificity),]


## Memory T cells
Th1.M.signature <- data.frame(gene=rownames(M.specs$specificities[M.specs$p.adj[,"Th1"] < fdr.Thres & M.specs$specificities[,"Th1"] > sp.Thres,]),
                              specificity=M.specs$specificities[M.specs$p.adj[,"Th1"] < fdr.Thres & M.specs$specificities[,"Th1"] > sp.Thres,]$Th1,
                              p.value=M.specs$p.val[M.specs$p.adj[,"Th1"] < fdr.Thres & M.specs$specificities[,"Th1"] > sp.Thres,]$Th1,
                              p.adj=M.specs$p.adj[M.specs$p.adj[,"Th1"] < fdr.Thres & M.specs$specificities[,"Th1"] > sp.Thres,]$Th1)

Th17.iTreg.M.signature <- data.frame(gene=rownames(M.specs$specificities[M.specs$p.adj[,"Th17.iTreg"] < fdr.Thres & M.specs$specificities[,"Th17.iTreg"] > sp.Thres,]),
                               specificity=M.specs$specificities[M.specs$p.adj[,"Th17.iTreg"] < fdr.Thres & M.specs$specificities[,"Th17.iTreg"] > sp.Thres,]$Th17.iTreg,
                               p.value=M.specs$p.val[M.specs$p.adj[,"Th17.iTreg"] < fdr.Thres & M.specs$specificities[,"Th17.iTreg"] > sp.Thres,]$Th17.iTreg,
                               p.adj=M.specs$p.adj[M.specs$p.adj[,"Th17.iTreg"] < fdr.Thres & M.specs$specificities[,"Th17.iTreg"] > sp.Thres,]$Th17.iTreg)

IFNB.M.signature <- data.frame(gene=rownames(M.specs$specificities[M.specs$p.adj[,"IFNB"] < fdr.Thres & M.specs$specificities[,"IFNB"] > sp.Thres,]),
                               specificity=M.specs$specificities[M.specs$p.adj[,"IFNB"] < fdr.Thres & M.specs$specificities[,"IFNB"] > sp.Thres,]$IFNB,
                               p.value=M.specs$p.val[M.specs$p.adj[,"IFNB"] < fdr.Thres & M.specs$specificities[,"IFNB"] > sp.Thres,]$IFNB,
                               p.adj=M.specs$p.adj[M.specs$p.adj[,"IFNB"] < fdr.Thres & M.specs$specificities[,"IFNB"] > sp.Thres,]$IFNB)


Th1.M.signature <- Th1.M.signature[order(-Th1.M.signature$specificity),]
Th17.iTreg.M.signature <- Th17.iTreg.M.signature[order(-Th17.iTreg.M.signature$specificity),]
IFNB.M.signature <- IFNB.M.signature[order(-IFNB.M.signature$specificity),]

# SAVING RESULTS AS CSV FILES
write.csv(TH1.N.signature, "../Results/proteogenomic_signatures/N_Th1_signature.csv",quote=F,row.names=F)
write.csv(TH2.N.signature, "../Results/proteogenomic_signatures/N_Th2_signature.csv",quote=F,row.names=F)
write.csv(TH17.N.signature, "../Results/proteogenomic_signatures/N_Th17_signature.csv",quote=F,row.names=F)
write.csv(ITREG.N.signature, "../Results/proteogenomic_signatures/N_iTreg_signature.csv",quote=F,row.names=F)
write.csv(IFNB.N.signature, "../Results/proteogenomic_signatures/N_IFNB_signature.csv",quote=F,row.names=F)

write.csv(TH1.M.signature, "../Results/proteogenomic_signatures/M_Th1_signature.csv",quote=F,row.names=F)
write.csv(TH17.ITREG.M.signature, "../Results/proteogenomic_signatures/M_Th17_iTreg_signature.csv",quote=F,row.names=F)
write.csv(IFNB.M.signature, "../Results/proteogenomic_signatures/M_IFNB_signature.csv",quote=F,row.names=F)
