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
bioPkgTest("GeneOverlap")
bioPkgTest("edgeR")
pkgTest("VennDiagram")
pkgTest("gplots")

########################################################
#                   Gene Overlap                       #
########################################################

# Load non-normalized single cell experiment object
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

# RPKM normalize the data based on the total number of genes
x_rpkm <- rpkm(assay(sceset, "counts"), 28692)
rpkmsum <- sum(x_rpkm, na.rm=F)
tpm.values <- x_rpkm/rpkmsum * 10^6
assay(sceset, "RPKM") <- tpm.values

saveRDS(sceset, "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_rpkm.rds")

# Put AD and WT in different datasets
APP <- filter(sceset, Condition == "APP")
WT <- filter(sceset, Condition == "WT")

# Create dataframe for all nuclei in the AD dataframe
APP.nuclei <- filter(APP, Type == "Nuclei")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(APP.nuclei) > 0) > 0
APP.nuclei <- APP.nuclei[keep_feature, ]

# Create dataframe for all cells in the AD dataframe
APP.cells <- filter(APP, Type == "Cells")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(APP.cells) > 0) > 0
APP.cells <- APP.cells[keep_feature, ]

# Create dataframe for all nuclei in the WT dataframe
WT.nuclei <- filter(WT, Type == "Nuclei")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(WT.nuclei) > 0) > 0
WT.nuclei <- WT.nuclei[keep_feature, ]

# Create dataframe for all cells in the WT dataframe
WT.cells <- filter(WT, Type == "Cells")
# remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(WT.cells) > 0) > 0
WT.cells <- WT.cells[keep_feature, ]

# Create gene overlap object for AD
# Genome size is the size of the total amount of genes in the sceset
overlapAPP.obj <- newGeneOverlap(rownames(APP.nuclei), rownames(APP.cells), genome.size = 28692)

# Perform Fisher's exact test and calculate the Jaccard index to see if overlap is significant
overlapAPP.obj <- testGeneOverlap(overlapAPP.obj)

# Define the overlapping genes in AD
overlapAPP <- overlapAPP.obj@intersection

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
APP.cells.NO <- length(rownames(APP.cells)) - length(overlapAPP)

# Calculate the number of genes that are unique in the AD nuclei
APP.nuclei.NO <- length(rownames(APP.nuclei)) - length(overlapAPP)

####################################
# 25th percentile genes overlap    #
####################################


# Calculate the 25th percentile for AD nuclei
APP.nuclei.percentile <- assay(APP.nuclei, "RPKM")[rowSums(assay(APP.nuclei, "RPKM") >0) > quantile(rowSums(assay(APP.nuclei, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for AD cells
APP.cells.percentile <- assay(APP.cells, "RPKM")[rowSums(assay(APP.cells, "RPKM") >0) > quantile(rowSums(assay(APP.cells, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for WT nuclei
WT.nuclei.percentile <- assay(WT.nuclei, "RPKM")[rowSums(assay(WT.nuclei, "RPKM") >0) > quantile(rowSums(assay(WT.nuclei, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for WT cells
WT.cells.percentile <- assay(WT.cells, "RPKM")[rowSums(assay(WT.cells, "RPKM") >0) > quantile(rowSums(assay(WT.cells, "RPKM")) , probs=0.75),]

# Put all genenames from the percentiles in a list
geneList <- list(rownames(APP.nuclei.percentile), rownames(APP.cells.percentile), rownames(WT.nuclei.percentile), rownames(WT.cells.percentile))

# Rename the vectors
names(geneList) <- c("APP.nuclei", "APP.cells", "WT.nuclei", "WT.cells")

# Create a venn diagram with the genenames of each group
venn.plot <- venn.diagram(geneList , NULL, cex = 2, category.names=c("APP.nuclei", "APP.cells", "WT.nuclei", "WT.cells"), main="Overlap of the 25th percentile genes in each group", main.cex = 1.5)

# Plot the venn diagram
grid.draw(venn.plot)

# Get the list of genes present in each Venn compartment with the use of the gplot package
genes <- venn(geneList, show.plot=FALSE)

# Get the intersections
inters <- attr(genes, "intersections")

# Unique genes of each cell type and condition
WT.nuclei.df <- sceset[rownames(sceset) %in% inters$WT.nuclei, ]

# Get the gene names of the unique genes
rownames(WT.nuclei.df)

WT.cells.df <- sceset[rownames(sceset) %in% inters$WT.cells, ]

rownames(WT.cells.df)

AD.nuclei.df <- sceset[rownames(sceset) %in% inters$AD.nuclei, ]

rownames(AD.nuclei.df)

AD.cells.df <- sceset[rownames(sceset) %in% inters$AD.cells, ]

rownames(AD.cells.df)
