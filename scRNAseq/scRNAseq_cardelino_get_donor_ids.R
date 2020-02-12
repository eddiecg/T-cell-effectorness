### CODE TO GENERATE DONOR IDs FOR EACH SINGLE-CELL BASED ON NATURAL GENETIC VARIATION ###
## Author: Yuanhua Huang (yuanhua@ebi.ac.uk)
##
## This script (now deprecated and replaced by the Python implementation 'Vireo') was kindly shared to us by Dr Oliver Stegle's group
## It estimates posterior probabilities of cells belonging to a given donor based on the results from cellSNP

# LOADING LIBRARIES
library(vcfR)
library(cardelino)

# FETCHING INPUT ARGUMENTS
args <- commandArgs(TRUE)
vcf_file <- args[1]
n_donors <- as.integer(args[2])

if (length(args) >= 3){
  out_file <- args[3]
} else{
  out_file <- "donorID_res.rds"
}
if (length(args) >= 4){
  check_doublet <- as.logical(args[4])
} else{
  check_doublet <- TRUE
}


# DEFINING FUNCTIONS
donor_id_out <- function(out, A, D, p_threshold=0.95, check_doublet=TRUE,
                         n_vars_threshold=10){
  ## output data
  out$A <- A
  out$D <- D
  # out$GT <- in_data$GT_cells #out has GT output
  
  ## assign data frame
  n_vars <- matrixStats::colSums2(!is.na(out$D))
  assigned <- assign_cells_to_clones(out$prob, threshold = p_threshold)
  colnames(assigned) <- c("cell", "donor_id", "prob_max")
  if (check_doublet) {
    assigned$prob_doublet <- matrixStats::rowSums2(out$prob_doublet)
    assigned$donor_id[assigned$prob_doublet > 0.05] <- "doublet"
  } else {
    assigned$prob_doublet <- NA
  }
  assigned$n_vars <- n_vars
  assigned$donor_id[n_vars < n_vars_threshold] <- "unassigned"
  assigned$donor_id[assigned$prob_doublet < 0.5 & rowSums(out$prob) <0.95] <- "unassigned"
  
  out$assigned <- assigned
  out
}

# ESTIMATING DONOR IDs
vcf_temp <- read.vcfR(vcf_file)
dp_full <- extract.gt(vcf_temp, element = "DP", as.numeric = TRUE)
ad_full <- extract.gt(vcf_temp, element = "AD", as.numeric = TRUE)

dp_sum <- extract.info(vcf_temp, element = "DP", as.numeric=TRUE)
oth_sum <- extract.info(vcf_temp, element = "OTH", as.numeric=TRUE)
sum(oth_sum / dp_sum > 0.2)

idx <- ((oth_sum/dp_sum < 0.2) & rowSums(dp_full, na.rm = TRUE) >= 20 & 
          (rowSums(ad_full, na.rm = TRUE)/rowSums(dp_full, na.rm = TRUE)>0.1))
mean(idx)
sum(idx)
A = ad_full[idx, ]
D = dp_full[idx, ]

don_ids <- donor_id_EM(A=A, D=D, K = n_donors, n_init=10, check_doublet=check_doublet)
ids <- donor_id_out(don_ids, A = A, D = D, check_doublet=check_doublet)

# EXPORTING RESULTS AS .RDS FILES
saveRDS(ids, out_file)