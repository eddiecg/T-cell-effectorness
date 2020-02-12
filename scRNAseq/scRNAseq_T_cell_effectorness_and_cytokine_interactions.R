### MODELLING GENE EXPRESSION AS A FUNCTION OF EFFECTORNESS AND CYTOKINES ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code models gene expression as a function of T cell effectorness and cytokines
# The expression values of a set of genes of interest is fitted into a linear regression model with the following design:
#
# gene_expression ~ T_cell_effectorness + cytokine_condition + T_cell_effectorness:cytokine_condition
#
# This model equation is used to estimate regression coefficients for each term (including the interaction term). 
# P values are then estimated using ANOVA and corrected for multiple testing using FDR

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
## Loading a Seurat object containing the expression matrix, annotations, and effectorness values for all cells in the study
## Please refer to "scRNAseq_T_cell_effectornes_analysis.R" to see how this object was generated
T.cells <- readRDS("allCells_Seurat_with_effectorness.rds")

# EXTRACTING STIMULATED CELLS
stim.cells <- SubsetData(T.cstim.cellsells, cells.use = T.cells@cell.names[T.cells@meta.data$cytokine.condition!="UNS"])

# IDENTIFYING GENES OF INTEREST
## Identifying highly variable genes
stim.cells <- ScaleData(stim.cells, vars.to.regress = c("S.Score","G2M.Score","donor.id","percent.mito","nUMI"))
stim.cells <- FindVariableGenes(stim.cells,x.low.cutoff=0.1,x.high.cutoff=5,y.cutoff=2,do.plot=F)

## Extracting genes of interest (highly variable genes)
genelist <- stim.cells@var.genes


# MODELLING GENE EXPRESSION AS A FUNCTION OF EFFECTORNESS
## Building a gene expression linear model, with effectorness and cytokine condition as independent variables and allowing for an interaction between them
## Note that effectorness is modeled as a continuous variable, while cytokine condition is a categorical variable
res <- sapply(genelist, FUN=function(genename){
  gene <- data.frame(exp=cbind(stim.cells@data[genename,]),eff=stim.cells@meta.data$effectorness, cyt=stim.cells@meta.data$cytokine.condition)
  gene <- gene[gene$exp > 0,]
  g <- lm(exp ~ eff*cyt, data=gene)
  res <- cbind(coeff=g$coefficients, pval=summary(g)$coefficients[,4])
  return(res)
})

## Converting regression results to the correct format
res <- data.frame(t(res))
colnames(res) <- c(paste("coeff",c("Int","eff","cyt.Th2","cyt.Th17","cyt.iTreg","eff:cyt.Th2","eff:cyt.Th17","eff:cyt.iTreg"),sep="_"),
                   paste("pval",c("Int","eff","cyt.Th2","cyt.Th17","cyt.iTreg","eff:cyt.Th2","eff:cyt.Th17","eff:cyt.iTreg"),sep="_"))

## Adjusting P values for multiple testing using FDR
padj <- apply(res[,grep("pval",colnames(res))],MARGIN=2,FUN=function(p){
  p.adjust(p, method="fdr")
})
colnames(padj) <- gsub("pval","padj",colnames(padj))
res <- cbind(res,padj)

## Performing ANOVA testing
ps <- sapply(genelist, FUN=function(genename){
  gene <- data.frame(exp=cbind(stim.cells@data[genename,]),eff=stim.cells@meta.data$effectorness, cyt=stim.cells@meta.data$cytokine.condition)
  gene <- gene[gene$exp > 0,]
  g <- lm(exp ~ eff*cyt, data=gene)
  res <- anova(g)$`Pr(>F)`[1:3]
  return(res)
})
ps <- t(ps)

ps <- apply(ps,MARGIN=2,FUN=function(p){
  padj <- p.adjust(p, method="fdr")
  return(padj)
})
ps[ps[,3] > 0.01,]

## Saving results
write.csv(res, "effectorness_and_cytokine_linear_regression_coefficients.csv",quote=F)
write.csv(ps, "effectorness_anova_results.csv",quote=F)