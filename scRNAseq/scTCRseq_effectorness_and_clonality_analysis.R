### ANALYSIS OF THE RELATIONSHIP BETWEEN T CELL CLONALITY AND T CELL EFFECTORNESS IN PUBLICLY AVAILABLE DATA  ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code leverages a public data set containing paired TCR-seq data and scRNA-seq data for T cells from blood and tumours in colorectal cancer
# First, the data is explored and cells are ordered into effectorness trajectories using monocle 2 (see "scRNAseq_T_cell_effectorness_analysis.R" for more information)
# Next, the clonotypes of T cells at different ranges of effectorness are analysed to assess their diversity, as well as Shannon entropy and other metrics

# LOADING LIBRARIES
library(ggplot2)
library(Seurat)
library(monocle)
library(zoo)

# LOADING DATA
## Loading public TCR-seq + scRNA-seq. These data were obtained from https://doi.org/10.1038/s41586-018-0694-x
cellsPath = gzfile("./GSE108989_CRC.TCell.S11138.count.txt.gz","rt")  
cellsCounts = read.csv(cellsPath,header=T,sep="\t")

cellsMetadata = read.csv("GSE108989_CRC.TCell.S11138.TCRtyping.csv", header=T)

# FORMATING DATA
## Formatting counts table
cellsCounts = cellsCounts[,-1]
cellsCounts = cellsCounts[!is.na(cellsCounts$symbol),]
rownames(cellsCounts) = cellsCounts$symbol
cellsCounts = cellsCounts[,-1]
colnames(cellsCounts) = gsub("\\.","-",colnames(cellsCounts))

## Subsetting counts table
cellNames = intersect(cellsMetadata$Cell.name, colnames(cellsCounts))
cellsCounts = cellsCounts[,cellNames]

## Formatting metadata table
rownames(cellsMetadata) = cellsMetadata$Cell.name
cellsMetadata = cellsMetadata[cellNames,]
colnames(cellsMetadata) = gsub("-",".",gsub("--",".",gsub("\\.","-",gsub("\\.$","",colnames(cellsMetadata)))))
colnames(cellsMetadata) = gsub("statusd","status",gsub("Clusterc","Cluster",gsub("geneb","gene",gsub("typea","type",colnames(cellsMetadata)))))

## Creating Seurat object
cells <- CreateSeuratObject(raw.data = cellsCounts, min.cells = 10, min.genes = 200, project = "CRC_TCR_data")

## Adding cell metadata
cells@meta.data = cbind(cells@meta.data, cellsMetadata)
colnames(cells@meta.data)[2] = "nCounts"

## Calculating MT content
estimateMitochondrial = function(dat){
  mito.genes = grep(pattern = "^MT", x=rownames(dat@data), value=TRUE)
  percent.mito = colSums(as.matrix(dat@raw.data[mito.genes,]))/colSums(as.matrix(dat@raw.data))
  dat = AddMetaData(dat, metadata = percent.mito, col.name = "percent.mito")
  return(dat)
}
cells = estimateMitochondrial(cells)

# PRE-PROCESSING DATA
## Keeping only CD4 T cells
cellsCD4 = SubsetData(cells, cells.use = cells@cell.names[cells@meta.data$Subtype == "CD4"])

## Adding source information
cellsCD4@meta.data$Tissue = NA
cellsCD4@meta.data$Tissue[grep("NT",cellsCD4@meta.data$Cell.type)] = "Adjacent.normal"
cellsCD4@meta.data$Tissue[grep("PT",cellsCD4@meta.data$Cell.type)] = "Tumour"
cellsCD4@meta.data$Tissue[grep("TT",cellsCD4@meta.data$Cell.type)] = "Peripheral.blood"

cellsCD4@meta.data$Tissue[grep("NP7",cellsCD4@meta.data$Cell.type)] = "Adjacent.normal"
cellsCD4@meta.data$Tissue[grep("TP7",cellsCD4@meta.data$Cell.type)] = "Tumour"
cellsCD4@meta.data$Tissue[grep("PP7",cellsCD4@meta.data$Cell.type)] = "Peripheral.blood"

# Keeping cells from peripheral blood only
cellsCD4blood = SubsetData(cellsCD4, cells.use = cellsCD4@cell.names[cellsCD4@meta.data$Tissue == "Peripheral.blood"])

# Removing Tregs (CD4+ CD25high)
cellsCD4blood = SubsetData(cellsCD4blood, cells.use = cellsCD4@cell.names[grep("TR",cellsCD4@meta.data$Cell.type, invert = T)])

# Normalizing data
cellsCD4blood <- NormalizeData(cellsCD4blood, normalization.method = "LogNormalize")

## Scoring cells based on their expression of cell cycle genes
cc.genes = readLines("regev_lab_cell_cycle_genes.txt")
S.genes = cc.genes[1:43]
G2M.genes = cc.genes[44:97]

cellsCD4blood = CellCycleScoring(cellsCD4blood, s.genes=S.genes, g2m.genes=G2M.genes, set.ident=TRUE)

# EXPLORATORY DATA ANALYSIS
## Reducing dimensions
cellsCD4blood = ScaleData(cellsCD4blood, vars.to.regress = c("nCounts","percent.mito","Patient", "Tissue"))
cellsCD4blood = FindVariableGenes(cellsCD4blood, x.low.cutoff = 0.1, x.high.cutoff = 4, y.cutoff = 2, do.plot = T)
cellsCD4blood = RunPCA(cellsCD4blood, pc.genes = cellsCD4blood@var.genes, do.print = FALSE, pcs.compute = 15)
cellsCD4blood = RunUMAP(cellsCD4blood, dim.use=1:15) #Run from a venv with umap!

## Visualising cells
DimPlot(cellsCD4blood, reduction.use = "umap", pt.size = 2, group.by = "Tissue")
DimPlot(cellsCD4blood, reduction.use = "umap", pt.size = 2, group.by = "Cluster")
FeaturePlot(cellsCD4blood, reduction.use = "umap", pt.size = 1, features.plot = c("SELL","CCR7","IL7R","KLRB1"), cols.use = c("lightgrey","darkblue"))
FeaturePlot(cellsCD4blood, reduction.use = "umap", pt.size = 1, features.plot = c("GZMA","GZMB","CCL4","CCL3"), cols.use = c("lightgrey","darkblue"))
FeaturePlot(cellsCD4blood, reduction.use = "umap", pt.size = 1, features.plot = c("FOXP3","CTLA4","IL2RA","IL7R"), cols.use = c("lightgrey","darkblue"))


# ORDERING CELLS IN A PSEUDOTIME TRAJECTORY
## Converting to Monocle format
cellsCD4bloodMonocle <- importCDS(cellsCD4blood)

## Processing data with Monocle pipeline
cellsCD4bloodMonocle <- estimateSizeFactors(cellsCD4bloodMonocle)
fData(cellsCD4bloodMonocle)$num_cells_expressed <- rowSums(as.matrix(exprs(cellsCD4bloodMonocle) > 0))

ordering_genes <- cellsCD4blood@var.genes
ordering_genes <- setdiff(ordering_genes, unique(unlist(cc.genes)))
cellsCD4bloodMonocle <- setOrderingFilter(cellsCD4bloodMonocle, ordering_genes)

cellsCD4bloodMonocle <- reduceDimension(cellsCD4bloodMonocle, max_components = 2,method = 'DDRTree')

## Ordering cells along trajectory
cellsCD4bloodMonocle <- orderCells(cellsCD4bloodMonocle)

# Inverting trajectory (based on visualisation)
cellsCD4bloodMonocle <- orderCells(cellsCD4bloodMonocle, root_state = 3)

## Identifying genes driving pseudotime
diff_test_res <- differentialGeneTest(cellsCD4bloodMonocle[ordering_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))

# CONVERTING PSEUDOTIME INTO EFFECTORNESS
## Definig effectorness as the pseudotime values of each cel scaled to the [0,1] range. 
## This make values comparable with the data set in our study
cellsCD4blood@meta.data$Pseudotime = cellsCD4bloodMonocle$Pseudotime
cellsCD4blood@meta.data$Effectorness = (pData(cellsCD4bloodMonocle)$Pseudotime - min(pData(cellsCD4bloodMonocle)$Pseudotime))/(max(pData(cellsCD4bloodMonocle)$Pseudotime)-min(pData(cellsCD4bloodMonocle)$Pseudotime))
cellsCD4blood@dr$traj = data.frame(t(cellsCD4bloodMonocle@reducedDimS))

cellsCD4blood@meta.data = data.frame(cellsCD4blood@meta.data)

# Visualizing trajectory
cellsDataFrame = data.frame(cbind(cellsCD4blood@dr$umap@cell.embeddings, cellsCD4blood@dr$traj, cellsCD4blood@meta.data))
cellsDataFrame$Cluster = droplevels(cellsDataFrame$Cluster)

ggplot(cellsDataFrame, aes(x=X1, y=X2)) + 
  geom_point(aes(color=Effectorness),size=2) + 
  scale_color_gradient(low = "#fed976",high = "#cc4c02") + 
  theme_classic() + 
  theme(legend.position = "bottom")

# CLONALITY ANALYSIS
## Estimating average clonality, diversity (i.e. Shannon entropy), and number of clones per cluster
cellsDataFrame$Cluster = factor(cellsDataFrame$Cluster, levels=c("CD4_C01-CCR7","CD4_C02-ANXA1","CD4_C03-GNLY","CD4_C10-FOXP3"))

numberOfClones = sapply(levels(cellsDataFrame$Cluster), FUN = function(cellType){
  clones = length(unique(cellsDataFrame[cellsDataFrame$Cluster==cellType,"Clone.ID"]))
  return(clones)
})

frequencyOfClones = sapply(levels(cellsDataFrame$Cluster), FUN = function(cellType){
  clones = length(unique(cellsDataFrame[cellsDataFrame$Cluster==cellType,"Clone.ID"]))/sum(cellsDataFrame$Cluster==cellType)
  return(clones)
})

shannonEntropy = function(cells){
  cells = cells[cells$Clone.ID!="",]
  cells$Clone.ID = droplevels(cells$Clone.ID)
  cloneFrequency = table(cells$Clone.ID)
  cloneRelativeFrequency = cloneFrequency/nrow(cells)
  entropy = -sum(cloneRelativeFrequency*log2(cloneRelativeFrequency))/log2(length(cloneRelativeFrequency))
  return(entropy)
}

entropy = sapply(levels(cellsDataFrame$Cluster),FUN=function(cellType){
  selectedCells = cellsDataFrame[cellsDataFrame$Cluster==cellType,]
  shannonEntropy(selectedCells)
})

## Estimating clonality for different T cell effectorness ranges
sortedData = cellsDataFrame[,c("Effectorness","Clone.status","Clone.ID","Frequency")]
sortedData = sortedData[sortedData$Clone.ID!="",]
sortedData = sortedData[order(sortedData$Effectorness),]

### Using a sliding window to estimate T cell clonality at increasing values of T cell effectorness
windowWidth=100
cellClonality = data.frame(cbind(effectorness = rollapply(data=sortedData$Effectorness, width=windowWidth, FUN=mean), 
                                 numberOfClones = rollapply(data=as.character(sortedData$Clone.ID), width=windowWidth, FUN=function(z){
                                   length(unique(z))
                                   }),
                                 clonality = rollapply(data=sortedData$Clone.status, width=windowWidth, FUN=function(z){
                                  sum(z=="Clonal")/windowWidth
                                  }),
                                 medianCloneFrequency = rollapply(data=sortedData$Frequency, width=windowWidth, FUN=function(z){
                                   median(z, na.rm = T)
                                 }),
                                 shannonEntropy = rollapply(data=as.character(sortedData$Clone.ID), width=windowWidth, FUN=function(z){
                                   cloneFrequency = table(z)
                                   cloneRelativeFrequency = cloneFrequency/windowWidth
                                   entropy = -sum(cloneRelativeFrequency*log2(cloneRelativeFrequency))/log2(length(cloneRelativeFrequency))
                                   return(entropy)
                                   })
                                 ))

## Visualising the relationship between clonality metrics and effectorness
ggplot(cellClonality, aes(x=effectorness, y=numberOfClones)) + geom_point(size=1) + theme_bw() + theme(panel.grid = element_blank())
ggplot(cellClonality, aes(x=effectorness, y=clonality)) + geom_point(size=1) + theme_bw() + theme(panel.grid = element_blank())
ggplot(cellClonality, aes(x=effectorness, y=shannonEntropy)) + geom_point(size=1) + theme_bw() + theme(panel.grid = element_blank())

