### USING THE fCI METHOD TO PERFORM DIFFERENIAL EXPRESSION ANALYSIS OF BULK RNA-SEQ AND PROTEOMICS DATA JOINTLY ###
## Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)

# LOADING LIBRARIES
library(fCI)
library(reshape2)
library(ggplot2)
library(RColorBrewer)


# PERFORMING fCI TESTING
## Defining a function which performs fCI in a pair of RNA-seq and proteomics samples
## This function loads a table with joint RNA-seq and proteomics data for each cell type and cytokine condition and performs differential testing
## Refer to "integration_of_RNAseq_and_proteomics_data.R" to see how these joint tables were generated
run_fCI <- function(datafile, controls, cases){
  fci=new("NPCI")
  targets=find.fci.targets(.Object = fci,
                           wt.indexes = controls,
                           df.indexes = cases,
                           data.file = datafile)
  Diff.Expr.Genes=show.targets(targets)
  Diff.Expr.Genes$Log2_FC <- as.numeric(as.character(Diff.Expr.Genes$Log2_FC))
  Diff.Expr.Genes <- Diff.Expr.Genes[order(-Diff.Expr.Genes$Log2_FC),]
  return(Diff.Expr.Genes)
}


## Running fCI on naive T cells
results_naive_Th0 <- run_fCI(datafile="naive_5d_Th0_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_naive_Th1 <- run_fCI(datafile="naive_5d_Th1_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_naive_Th2 <- run_fCI(datafile="naive_5d_Th2_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_naive_Th17 <- run_fCI(datafile="naive_5d_Th17_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_naive_iTreg <- run_fCI(datafile="naive_5d_iTreg_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_naive_IFNB <- run_fCI(datafile="naive_5d_IFNB_proteogenomic_data.txt",controls=list(1:3, 6:8),cases=list(4:5, 9:11))

# Running fCI on memory T cells
results_memory_Th0 <- run_fCI(datafile="memory_5d_Th0_proteogenomic_data.txt",controls=list(1:3, 7:8),cases=list(4:6, 9:11))
results_memory_Th1 <- run_fCI(datafile="memory_5d_Th1_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_memory_Th2 <- run_fCI(datafile="memory_5d_Th2_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_memory_Th17 <- run_fCI(datafile="memory_5d_Th17_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_memory_iTreg <- run_fCI(datafile="memory_5d_iTreg_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:12))
results_memory_IFNB <- run_fCI(datafile="memory_5d_IFNB_proteogenomic_data.txt",controls=list(1:3, 7:9),cases=list(4:6, 10:11))


# SAVING fCI RESULTS
## Naive T cells
write.csv(results_naive_Th0, file="naive_Th0_vs_resting_fCI_results.csv", quote = F, row.names = F)
write.csv(results_naive_Th2, file="naive_Th1_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_naive_Th2, file="naive_Th2_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_naive_Th17, file="naive_Th17_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_naive_iTreg, file="naive_iTreg_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_naive_IFNB, file="/naive_IFNB_vs_Th0_fCI_results.csv", quote = F, row.names = F)

## Memory T cells
write.csv(results_memory_Th0, file="memory_Th0_vs_resting_fCI_results.csv", quote = F, row.names = F)
write.csv(results_memory_Th1, file="memory_Th1_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_memory_Th2, file="memory_Th2_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_memory_Th17, file="memory_Th17_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_memory_iTreg, file="memory_iTreg_vs_Th0_fCI_results.csv", quote = F, row.names = F)
write.csv(results_memory_IFNB, file="memory_IFNB_vs_Th0_fCI_results.csv", quote = F, row.names = F)
