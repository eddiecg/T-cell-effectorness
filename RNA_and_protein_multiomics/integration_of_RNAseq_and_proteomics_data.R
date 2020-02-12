### MERGING RNA-SEQ COUNTS AND PROTEIN ABUNDANCES INTO A SINGLE DATA SET ###
## Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)

## LOADING LIBRARIES
library(DESeq2)

## LOADING DATA
### Loading protein abundances
proteomics_data <- readRDS("proteinAbundances_summarizedExperiment.rds")

### Loading RNA counts
rna_data <- readRDS("../RNAseq/rawCounts_bulkRNAseq_DESeq2.rds")

## PRE-PROCESSING DATA
### Separating samples by cell type
proteomics_data_naive <- proteomics_data[,proteomics_data$cell_type=="CD4_naive"]
proteomics_data_memory <- proteomics_data[,proteomics_data$cell_type=="CD4_memory"]

rna_data_naive <- rna_data[,rna_data$cell_type=="CD4_Naive"]
rna_data_memory <- rna_data[,rna_data$cell_type=="CD4_Memory"]

## NAIVE T CELL DATA
### Extracting relative abundances only (i.e. discarding metadata)
protein.names <- as.character(rowData(proteomics_data_naive)$Gene_name)
proteomics_data_naive <- assay(proteomics_data_naive)
rownames(proteomics_data_naive) <- protein.names

### Removing lowly expressed genes (i.e. at least 30 counts in at least 3 samples)
rna_data_naive <- rna_data_naive[rowSums(assay(rna_data_naive) > 30) > 3,]

### Extracting abundances for the 5-day time point only (i.e. same time point as in the proteomics experiment)
rna_data_naive <- rna_data_naive[,rna_data_naive$stimulation_time == "5d"]

### Extracting normalised RNA counts only (i.e. discarding metadata)
rna.names <- as.character(rowData(rna_data_naive)$Gene_symbol)
rna_data_naive <- counts(rna_data_naive, normalized=T)
rownames(rna_data_naive) <- rna.names

### Identifying genes detected in both assays (LC-MS/MS and RNA-seq)
commonGenes <- intersect(protein.names, rna.names)

### Subsetting RNA-seq and proteomics data
proteomics_data_naive <- proteomics_data_naive[commonGenes,]
rna_data_naive <- rna_data_naive[commonGenes,]

### Matching column names between data sets
colnames(proteomics_data_naive) <- paste("PROTEIN_", 
                                  gsub("resting","Resting",gsub("_CD4_naive","",colnames(proteomics_data_naive))), 
                                  sep="")
colnames(rna_data_naive) <- paste("RNA_D", 
                                  gsub("5d_","",gsub("_CD4_Naive","",colnames(rna_data_naive))), 
                                  sep="")

### Combining RNA-seq and proteomics data in a single data frame 
### This is done independently for each cytokine condition
naive_Th0_protein_and_RNA <- data.frame(rna_data_naive[,"RNA_D254_Resting"],
                                        rna_data_naive[,"RNA_D255_Resting"],
                                        rna_data_naive[,"RNA_D257_Resting"],
                                        rna_data_naive[,"RNA_D254_Th0"],
                                        rna_data_naive[,"RNA_D255_Th0"],
                                        rna_data_naive[,"RNA_D257_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D254_Resting"],
                                        proteomics_data_naive[,"PROTEIN_D255_Resting"],
                                        proteomics_data_naive[,"PROTEIN_D257_Resting"],
                                        proteomics_data_naive[,"PROTEIN_D254_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D255_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D257_Th0"])

naive_Th1_protein_and_RNA <- data.frame(rna_data_naive[,"RNA_D254_Th0"],
                                         rna_data_naive[,"RNA_D255_Th0"],
                                         rna_data_naive[,"RNA_D257_Th0"],
                                         rna_data_naive[,"RNA_D254_Th1"],
                                         rna_data_naive[,"RNA_D255_Th1"],
                                         rna_data_naive[,"RNA_D257_Th1"],
                                         proteomics_data_naive[,"PROTEIN_D254_Th0"],
                                         proteomics_data_naive[,"PROTEIN_D255_Th0"],
                                         proteomics_data_naive[,"PROTEIN_D257_Th0"],
                                         proteomics_data_naive[,"PROTEIN_D254_Th1"],
                                         proteomics_data_naive[,"PROTEIN_D255_Th1"],
                                         proteomics_data_naive[,"PROTEIN_D257_Th1"])

naive_Th2_protein_and_RNA <- data.frame(rna_data_naive[,"RNA_D254_Th0"],
                                        rna_data_naive[,"RNA_D255_Th0"],
                                        rna_data_naive[,"RNA_D257_Th0"],
                                        rna_data_naive[,"RNA_D254_Th2"],
                                        rna_data_naive[,"RNA_D255_Th2"],
                                        rna_data_naive[,"RNA_D257_Th2"],
                                        proteomics_data_naive[,"PROTEIN_D254_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D255_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D257_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D254_Th2"],
                                        proteomics_data_naive[,"PROTEIN_D255_Th2"],
                                        proteomics_data_naive[,"PROTEIN_D257_Th2"])
                  
naive_Th17_protein_and_RNA <- data.frame(rna_data_naive[,"RNA_D254_Th0"],
                                        rna_data_naive[,"RNA_D255_Th0"],
                                        rna_data_naive[,"RNA_D257_Th0"],
                                        rna_data_naive[,"RNA_D254_Th17"],
                                        rna_data_naive[,"RNA_D255_Th17"],
                                        rna_data_naive[,"RNA_D257_Th17"],
                                        proteomics_data_naive[,"PROTEIN_D254_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D255_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D257_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D254_Th17"],
                                        proteomics_data_naive[,"PROTEIN_D255_Th17"],
                                        proteomics_data_naive[,"PROTEIN_D257_Th17"])

naive_iTreg_protein_and_RNA <- data.frame(rna_data_naive[,"RNA_D254_Th0"],
                                        rna_data_naive[,"RNA_D255_Th0"],
                                        rna_data_naive[,"RNA_D257_Th0"],
                                        rna_data_naive[,"RNA_D254_iTreg"],
                                        rna_data_naive[,"RNA_D255_iTreg"],
                                        rna_data_naive[,"RNA_D257_iTreg"],
                                        proteomics_data_naive[,"PROTEIN_D254_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D255_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D257_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D254_iTreg"],
                                        proteomics_data_naive[,"PROTEIN_D255_iTreg"],
                                        proteomics_data_naive[,"PROTEIN_D257_iTreg"])

naive_IFNB_protein_and_RNA <- data.frame(rna_data_naive[,"RNA_D254_Th0"],
                                        rna_data_naive[,"RNA_D255_Th0"],
                                        rna_data_naive[,"RNA_D257_Th0"],
                                        rna_data_naive[,"RNA_D255_IFNB"],
                                        rna_data_naive[,"RNA_D257_IFNB"],
                                        proteomics_data_naive[,"PROTEIN_D254_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D255_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D257_Th0"],
                                        proteomics_data_naive[,"PROTEIN_D254_IFNB"],
                                        proteomics_data_naive[,"PROTEIN_D255_IFNB"],
                                        proteomics_data_naive[,"PROTEIN_D257_IFNB"])

### Trimming column names
colnames(naive_Th0_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_naive","",gsub("rna_data_naive","",colnames(naive_Th0_protein_and_RNA))))
colnames(naive_Th1_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_naive","",gsub("rna_data_naive","",colnames(naive_Th1_protein_and_RNA))))
colnames(naive_Th2_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_naive","",gsub("rna_data_naive","",colnames(naive_Th2_protein_and_RNA))))                  
colnames(naive_Th17_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_naive","",gsub("rna_data_naive","",colnames(naive_Th17_protein_and_RNA))))
colnames(naive_iTreg_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_naive","",gsub("rna_data_naive","",colnames(naive_iTreg_protein_and_RNA))))
colnames(naive_IFNB_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_naive","",gsub("rna_data_naive","",colnames(naive_IFNB_protein_and_RNA))))

## SAVING RESULTS
write.table(naive_Th0_protein_and_RNA, "naive_5d_Th0_proteogenomic_data.txt", sep="\t", quote=F)
write.table(naive_Th1_protein_and_RNA, "naive_5d_Th1_proteogenomic_data.txt", sep="\t", quote=F)
write.table(naive_Th2_protein_and_RNA, "naive_5d_Th2_proteogenomic_data.txt", sep="\t", quote=F)
write.table(naive_Th17_protein_and_RNA, "naive_5d_Th17_proteogenomic_data.txt", sep="\t", quote=F)
write.table(naive_iTreg_protein_and_RNA, "naive_5d_iTreg_proteogenomic_data.txt", sep="\t", quote=F)
write.table(naive_IFNB_protein_and_RNA, "naive_5d_IFNB_proteogenomic_data.txt", sep="\t", quote=F)


## MEMORY T CELL DATA
### Extracting relative abundances only (i.e. discarding metadata)
protein.names <- as.character(rowData(proteomics_data_memory)$Gene_name)
proteomics_data_memory <- assay(proteomics_data_memory)
rownames(proteomics_data_memory) <- protein.names

### Removing lowly expressed genes (i.e. at least 30 counts in at least 3 samples)
rna_data_memory <- rna_data_memory[rowSums(assay(rna_data_memory) > 30) > 3,]

### Extracting abundances for the 5-day time point only (i.e. same time point as in the proteomics experiment)
rna_data_memory <- rna_data_memory[,rna_data_memory$stimulation_time == "5d"]

### Extracting normalised RNA counts only (i.e. discarding metadata)
rna.names <- as.character(rowData(rna_data_memory)$Gene_symbol)
rna_data_memory <- counts(rna_data_memory, normalized=T)
rownames(rna_data_memory) <- rna.names

### Identifying genes detected in both assays (LC-MS/MS and RNA-seq)
commonGenes <- intersect(protein.names, rna.names)

### Subsetting RNA-seq and proteomics data
proteomics_data_memory <- proteomics_data_memory[commonGenes,]
rna_data_memory <- rna_data_memory[commonGenes,]

### Matching column names between data sets
colnames(proteomics_data_memory) <- paste("PROTEIN_", 
                                         gsub("resting","Resting",gsub("_CD4_memory","",colnames(proteomics_data_memory))), 
                                         sep="")
colnames(rna_data_memory) <- paste("RNA_D", 
                                  gsub("5d_","",gsub("_CD4_Memory","",colnames(rna_data_memory))), 
                                  sep="")

### Combining RNA-seq and proteomics data in a single data frame 
### This is done independently for each cytokine condition
memory_Th0_protein_and_RNA <- data.frame(rna_data_memory[,"RNA_D262_Resting"],
                                        rna_data_memory[,"RNA_D264_Resting"],
                                        rna_data_memory[,"RNA_D265_Resting"],
                                        rna_data_memory[,"RNA_D262_Th0"],
                                        rna_data_memory[,"RNA_D264_Th0"],
                                        rna_data_memory[,"RNA_D265_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D262_Resting"],
                                        proteomics_data_memory[,"PROTEIN_D264_Resting"],
                                        proteomics_data_memory[,"PROTEIN_D262_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D264_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D265_Th0"])

memory_Th1_protein_and_RNA <- data.frame(rna_data_memory[,"RNA_D262_Th0"],
                                        rna_data_memory[,"RNA_D264_Th0"],
                                        rna_data_memory[,"RNA_D265_Th0"],
                                        rna_data_memory[,"RNA_D262_Th1"],
                                        rna_data_memory[,"RNA_D264_Th1"],
                                        rna_data_memory[,"RNA_D265_Th1"],
                                        proteomics_data_memory[,"PROTEIN_D262_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D264_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D265_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D262_Th1"],
                                        proteomics_data_memory[,"PROTEIN_D264_Th1"],
                                        proteomics_data_memory[,"PROTEIN_D265_Th1"])

memory_Th2_protein_and_RNA <- data.frame(rna_data_memory[,"RNA_D262_Th0"],
                                        rna_data_memory[,"RNA_D264_Th0"],
                                        rna_data_memory[,"RNA_D265_Th0"],
                                        rna_data_memory[,"RNA_D262_Th2"],
                                        rna_data_memory[,"RNA_D264_Th2"],
                                        rna_data_memory[,"RNA_D265_Th2"],
                                        proteomics_data_memory[,"PROTEIN_D262_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D264_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D265_Th0"],
                                        proteomics_data_memory[,"PROTEIN_D262_Th2"],
                                        proteomics_data_memory[,"PROTEIN_D264_Th2"],
                                        proteomics_data_memory[,"PROTEIN_D265_Th2"])

memory_Th17_protein_and_RNA <- data.frame(rna_data_memory[,"RNA_D262_Th0"],
                                         rna_data_memory[,"RNA_D264_Th0"],
                                         rna_data_memory[,"RNA_D265_Th0"],
                                         rna_data_memory[,"RNA_D262_Th17"],
                                         rna_data_memory[,"RNA_D264_Th17"],
                                         rna_data_memory[,"RNA_D265_Th17"],
                                         proteomics_data_memory[,"PROTEIN_D262_Th0"],
                                         proteomics_data_memory[,"PROTEIN_D264_Th0"],
                                         proteomics_data_memory[,"PROTEIN_D265_Th0"],
                                         proteomics_data_memory[,"PROTEIN_D262_Th17"],
                                         proteomics_data_memory[,"PROTEIN_D264_Th17"],
                                         proteomics_data_memory[,"PROTEIN_D265_Th17"])

memory_iTreg_protein_and_RNA <- data.frame(rna_data_memory[,"RNA_D262_Th0"],
                                          rna_data_memory[,"RNA_D264_Th0"],
                                          rna_data_memory[,"RNA_D265_Th0"],
                                          rna_data_memory[,"RNA_D262_iTreg"],
                                          rna_data_memory[,"RNA_D264_iTreg"],
                                          rna_data_memory[,"RNA_D265_iTreg"],
                                          proteomics_data_memory[,"PROTEIN_D262_Th0"],
                                          proteomics_data_memory[,"PROTEIN_D264_Th0"],
                                          proteomics_data_memory[,"PROTEIN_D265_Th0"],
                                          proteomics_data_memory[,"PROTEIN_D262_iTreg"],
                                          proteomics_data_memory[,"PROTEIN_D264_iTreg"],
                                          proteomics_data_memory[,"PROTEIN_D265_iTreg"])

memory_IFNB_protein_and_RNA <- data.frame(rna_data_memory[,"RNA_D262_Th0"],
                                         rna_data_memory[,"RNA_D264_Th0"],
                                         rna_data_memory[,"RNA_D265_Th0"],
                                         rna_data_memory[,"RNA_D264_IFNB"],
                                         rna_data_memory[,"RNA_D265_IFNB"],
                                         proteomics_data_memory[,"PROTEIN_D262_Th0"],
                                         proteomics_data_memory[,"PROTEIN_D264_Th0"],
                                         proteomics_data_memory[,"PROTEIN_D265_Th0"],
                                         proteomics_data_memory[,"PROTEIN_D264_IFNB"],
                                         proteomics_data_memory[,"PROTEIN_D265_IFNB"])

### Trimming column names
colnames(memory_Th0_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_memory","",gsub("rna_data_memory","",colnames(memory_Th0_protein_and_RNA))))
colnames(memory_Th1_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_memory","",gsub("rna_data_memory","",colnames(memory_Th1_protein_and_RNA))))
colnames(memory_Th2_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_memory","",gsub("rna_data_memory","",colnames(memory_Th2_protein_and_RNA))))                  
colnames(memory_Th17_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_memory","",gsub("rna_data_memory","",colnames(memory_Th17_protein_and_RNA))))
colnames(memory_iTreg_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_memory","",gsub("rna_data_memory","",colnames(memory_iTreg_protein_and_RNA))))
colnames(memory_IFNB_protein_and_RNA) <- gsub("\\.","",gsub("proteomics_data_memory","",gsub("rna_data_memory","",colnames(memory_IFNB_protein_and_RNA))))

## SAVING RESULTS
write.table(memory_Th0_protein_and_RNA, "memory_5d_Th0_proteogenomic_data.txt", sep="\t", quote=F)
write.table(memory_Th1_protein_and_RNA, "memory_5d_Th1_proteogenomic_data.txt", sep="\t", quote=F)
write.table(memory_Th2_protein_and_RNA, "memory_5d_Th2_proteogenomic_data.txt", sep="\t", quote=F)
write.table(memory_Th17_protein_and_RNA, "memory_5d_Th17_proteogenomic_data.txt", sep="\t", quote=F)
write.table(memory_iTreg_protein_and_RNA, "memory_5d_iTreg_proteogenomic_data.txt", sep="\t", quote=F)
write.table(memory_IFNB_protein_and_RNA, "memory_5d_IFNB_proteogenomic_data.txt", sep="\t", quote=F)