### INTEGRATING RAW MASS SPEC ABUNANCES, PROTEIN ANNOTATIONS, AND SAMPLE ANNOTATIONS INTO A SINGLE OBJECT ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)

# LOADING LIBRARIES
library(SummarizedExperiment)
library(annotables)

# LOADING DATA
## Loading scaled abundances from Mass Spectrometry
MassSpec_data <- read.table("NCOMMS-19-7936188_MassSpec_scaled_abundances.txt", header = T, stringsAsFactors = F)

# CREATING SUMMARISED EXPERIMENT
### Creating a data frame with protein annotations
protein_annotations <- data.frame(MassSpec_data[,c("Protein_id","Gene_id","Gene_name")], row.names = MassSpec_data$Gene_name)

rownames(MassSpec_data) <- MassSpec_data$Gene_name
MassSpec_data <- MassSpec_data[,-c(1:3)]

### Creating a data frame with sample annotations
sample_ids <- colnames(MassSpec_data)
sample_annotations <- data.frame(row.names = sample_ids,
                                 donor_id = sapply(sample_ids, function(x){strsplit(x, split = "_")[[1]][1]}),
                                 cell_type = paste("CD4", 
                                                   sapply(sample_ids, function(x){strsplit(x, split = "_")[[1]][3]}), 
                                                   sep="_"),
                                 cytokine_condition = sapply(sample_ids, function(x){strsplit(x, split = "_")[[1]][4]}),
                                 stringsAsFactors = T)
sample_annotations$activation_status <- ifelse(sample_annotations$cytokine_condition == "resting", "Resting", "Activated")

### Creating a variable with relevant metadata for the study
meta <- list(
  Study="Mapping cytokine induced gene expression changes in human CD4+ T cells",
  Experiment="Quantitative proteomics (LC-MS/MS) panel of cytokine induced T cell polarisations",
  Laboratory="Trynka Group, Wellcome Sanger Institute",
  Experimenter=c("Eddie Cano-Gamez",
                 "Blagoje Soskic",
                 "Deborah Plowman"),
  Description="To study cytokine-induced cell polarisation, we isolated human naive and memory CD4+ T cells in triplicate from peripheral blood of healthy individuals. Next, we polarised the cells with different cytokine combinations linked to autoimmunity and performed LC-MS/MS.",
  Methdology="LC-MS/MS with isobaric labelling",
  Characteristics="Data type: Normalised, scaled protein abundances",
  Date="September, 2019",
  URL="https://doi.org/10.1101/753731"
)

### Creating summarised experiment object
proteomics_data <- SummarizedExperiment(assays=list(counts=as.matrix(MassSpec_data)),
                                        colData=sample_annotations, 
                                        rowData=protein_annotations, 
                                        metadata=meta)


# SAVING SUMMARISED EXPERIMENT AS BINARY FILE
saveRDS(proteomics_data, file="proteinAbundances_summarizedExperiment.rds")
