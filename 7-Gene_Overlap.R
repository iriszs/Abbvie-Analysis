# Remove all environmental variables
rm(list=ls(all=TRUE))
# Restart R to make sure a 'clean' environment is created and all used memory is freed
.rs.restartR()

########################################################
#                      Sources                         #
########################################################

# Load script that contains the used functions in this script
source("/media/imgorter/BD_1T/Iris/Scripts/abbvie/workflow/packageTest.R")

########################################################
#                     Options                          #
########################################################

options(stringsAsFactors = FALSE)
# Pseudo-random number generation to insure results can be reproduced
set.seed(1234567)

########################################################
#                     Packages                         #
########################################################

# Load only the necessary packages for this script
# With the use of the pkgTest and bioPkgTest functions from the packageTest.r file
bioPkgTest("SingleCellExperiment")
bioPkgTest("clusterProfiler")
bioPkgTest("org.Mm.eg.db")
bioPkgTest("scater")
pkgTest("VennDiagram")
pkgTest("gplots")

########################################################
#                   Gene Overlap                       #
########################################################

# Load non-normalized single cell experiment object
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

# Put AD and WT in different datasets
AD <- filter(sceset, Condition == "AD")
WT <- filter(sceset, Condition == "WT")

# Remove the single cell experiment object
rm("sceset")
gc()

# Create dataframe for all nuclei in the AD dataframe
AD.nuclei <- filter(AD, Type == "Nuclei")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(AD.nuclei) > 0) > 0
AD.nuclei <- AD.nuclei[keep_feature, ]

# Create dataframe for all cells in the AD dataframe
AD.cells <- filter(AD, Type == "Microglia")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(AD.cells) > 0) > 0
AD.cells <- AD.cells[keep_feature, ]

# Create dataframe for all nuclei in the WT dataframe
WT.nuclei <- filter(WT, Type == "Nuclei")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(WT.nuclei) > 0) > 0
WT.nuclei <- WT.nuclei[keep_feature, ]

# Create dataframe for all cells in the WT dataframe
WT.cells <- filter(WT, Type == "Microglia")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(WT.cells) > 0) > 0
WT.cells <- WT.cells[keep_feature, ]

# Create gene overlap object for AD
# Genome size is the size of the total amount of genes in the sceset
overlapAD.obj <- newGeneOverlap(rownames(AD.nuclei), rownames(AD.cells), genome.size = 28692)

# Perform Fisher's exact test and calculate the Jaccard index to see if overlap is significant
overlapAD.obj <- testGeneOverlap(overlapAD.obj)

# Define the overlapping genes in AD
overlapAD <- overlapAD.obj@intersection

# Create gene overlap object for WT
# Genome size is the size of the total amount of genes in the sceset
overlapWT.obj <- newGeneOverlap(rownames(WT.nuclei), rownames(WT.cells), genome.size = 28692)

# Perform Fisher's exact test and calculate the Jaccard index to see if overlap is significant
overlapWT.obj <- testGeneOverlap(overlapWT.obj)

# Define the overlapping genes in WT
overlapWT <- overlapWT.obj@intersection

# Calculate the number of genes that are unique in the WT cells
WT.cells.NO <- length(rownames(WT.cells)) - length(overlapWT)

# Calculate the number of genes that are unique in the WT nuclei
WT.nuclei.NO <- length(rownames(WT.nuclei)) - length(overlapWT)

# Calculate the number of genes that are unique in the AD cells
AD.cells.NO <- length(rownames(AD.cells)) - length(overlapAD)

# Calculate the number of genes that are unique in the AD nuclei
AD.nuclei.NO <- length(rownames(AD.nuclei)) - length(overlapAD)

####################################
# 25th percentile genes overlap    #
####################################


# Calculate the 25th percentile for AD nuclei
AD.nuclei.percentile <- assay(AD.nuclei, "RPKM")[rowSums(assay(AD.nuclei, "RPKM") >0) > quantile(rowSums(assay(AD.nuclei, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for AD cells
AD.cells.percentile <- assay(AD.cells, "RPKM")[rowSums(assay(AD.cells, "RPKM") >0) > quantile(rowSums(assay(AD.cells, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for WT nuclei
WT.nuclei.percentile <- assay(WT.nuclei, "RPKM")[rowSums(assay(WT.nuclei, "RPKM") >0) > quantile(rowSums(assay(WT.nuclei, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for WT cells
WT.cells.percentile <- assay(WT.cells, "RPKM")[rowSums(assay(WT.cells, "RPKM") >0) > quantile(rowSums(assay(WT.cells, "RPKM")) , probs=0.75),]

# Put all genenames from the percentiles in a list
geneList <- list(rownames(AD.nuclei.percentile), rownames(AD.cells.percentile), rownames(WT.nuclei.percentile), rownames(WT.cells.percentile))

# Rename the vectors
names(geneList) <- c("AD.nuclei", "AD.cells", "WT.nuclei", "WT.cells")

# Create a venn diagram with the genenames of each group
venn.plot <- venn.diagram(geneList , NULL, cex = 2, category.names=c("AD.nuclei", "AD.cells", "WT.nuclei", "WT.cells"), main="Overlap of the 25th percentile genes in each group", main.cex = 1.5)

# Plot the venn diagram
grid.draw(venn.plot)

# Get the list of genes present in each Venn compartment with the use of the gplot package
genes <- venn(geneList, show.plot=FALSE)

# unique genes of each cell type and condition
WT.nuclei.df <- sceset[rownames(sceset) %in% inters$WT.nuclei, ]

WT.cells.df <- sceset[rownames(sceset) %in% inters$WT.cells, ]

AD.nuclei.df <- sceset[rownames(sceset) %in% inters$AD.nuclei, ]

AD.cells.df <- sceset[rownames(sceset) %in% inters$AD.cells, ]

