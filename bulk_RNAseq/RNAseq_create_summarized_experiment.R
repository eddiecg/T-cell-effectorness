### INTEGRATING RAW RNA-SEQ COUNTS, GENE ANNOTATIONS, AND SAMPLE ANNOTATIONS INTO A SINGLE OBJECT ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)

# LOADING LIBRARIES
library(annotables)
library(SummarizedExperiment)

# FORMATING DATA
## Loading counts table (obtained from featureCounts)
counts <- read.table("NCOMMS-19-7936188_bulk_RNAseq_raw_counts.txt", header=T, row.names=1)

## Loading sample metadata
sample.info <- read.table("NCOMMS-19-7936188_bulk_RNAseq_metadata.txt", header=T, row.names=1, stringsAsFactors = T)

# FORMATING DATA
## Converting all annotations to the mos suitable data types
sample.info$sequencing_batch <- factor(sample.info$sequencing_batch)
sample.info$cell_culture_batch <- factor(sample.info$cell_culture_batch)
sample.info$activation_status <- ifelse(sample.info$cytokine_condition == "Resting", "Resting", "Activated")

## Fetching gene annotations from Ensembl87/GRCh38
features.data <- data.frame(grch38)
features.data <- features.data[!duplicated(features.data[c("ensgene")]),]
rownames(features.data) <- features.data$ensgene
features.data <- features.data[,-1]

gene.info <- data.frame(Gene_id = rownames(counts),
                        Gene_symbol = features.data[rownames(counts),"symbol"],
                        Chr = features.data[rownames(counts),"chr"],
                        Start = features.data[rownames(counts),"start"],
                        End = features.data[rownames(counts),"end"],
                        Strand = features.data[rownames(counts),"strand"],
                        Biotype = features.data[rownames(counts),"biotype"])

## Converting gene annotations to the most suitable data types
gene.info$Gene_id <- as.character(gene.info$Gene_id)
gene.info$Start <- as.numeric(gene.info$Start)
gene.info$End <- as.numeric(gene.info$End)
gene.info$Strand <- factor(gene.info$Strand)

## Creating a variable with relevant metadata for the study
meta <- list(
  Study="Mapping cytokine induced gene expression changes in human CD4+ T cells",
  Experiment="RNA-seq panel of cytokine induced T cell polarisations",
  Laboratory="Trynka Group, Wellcome Sanger Institute",
  Experimenter=c("Eddie Cano-Gamez",
                 "Blagoje Soskic",
                 "Deborah Plowman"),
  Description="To study cytokine-induced cell polarisation, we isolated human naive and memory CD4+ T cells in triplicate from peripheral blood of healthy individuals. Next, we polarised the cells with different cytokine combinations linked to autoimmunity and performed RNA-sequencing.",
  Methdology=c("Library Prep: Illumina TruSeq (Poly-A capture method)", 
               "Sequencing Platform: Illumina HiSeq 2500",
               "Read Alignment: STAR (MAPQ > 20)",
               "Read Quantification: featureCounts",
               "Reference: Ensembl 87 (GRCh38/hg38)"),
  Characteristics="Data type: Raw counts",
  Date="September, 2019",
  URL="https://doi.org/10.1101/753731"
)


# CREATING SUMMARIZED EXPERIMENT OBJECT
RNA.experiment <- SummarizedExperiment(assays=list(counts=as.matrix(counts)),
                     colData=sample.info,
                     rowData=gene.info,
                     metadata=meta)

# SAVING RESULTS AS R DATA SET
saveRDS(RNA.experiment, "bulkRNAseq_summarizedExperiment.rds")
