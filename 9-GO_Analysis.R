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
#                     Packages                         #
########################################################

# Load the RDS with the ensembl symbols
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_ensembl_Normalized.rds")

# Seperate AD and WT in different datasets
AD <- filter(sceset, Condition == "AD")
WT <- filter(sceset, Condition == "WT")

# Remove the single cell experiment object
rm("sceset")
gc()

# Create dataframe for all nuclei in the AD dataframe
AD.nuclei <- filter(AD, Type == "Nuclei")

# Create dataframe for all cells in the AD dataframe
AD.cells <- filter(AD, Type == "Microglia")

# Create dataframe for all nuclei in the WT dataframe
WT.nuclei <- filter(WT, Type == "Nuclei")

# Create dataframe for all cells in the WT dataframe
WT.cells <- filter(WT, Type == "Microglia")

# Calculate the 25th percentile for AD nuclei
AD.nuclei.percentile <- assay(AD.nuclei, "RPKM")[rowSums(assay(AD.nuclei, "RPKM") >0) > quantile(rowSums(assay(AD.nuclei, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for AD cells
AD.cells.percentile <- assay(AD.cells, "RPKM")[rowSums(assay(AD.cells, "RPKM") >0) > quantile(rowSums(assay(AD.cells, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for WT nuclei
WT.nuclei.percentile <- assay(WT.nuclei, "RPKM")[rowSums(assay(WT.nuclei, "RPKM") >0) > quantile(rowSums(assay(WT.nuclei, "RPKM")) , probs=0.75),]

# Calculate the 25th percentile for WT cells
WT.cells.percentile <- assay(WT.cells, "RPKM")[rowSums(assay(WT.cells, "RPKM") >0) > quantile(rowSums(assay(WT.cells, "RPKM")) , probs=0.75),]

# Save the genenames that are in the 25th percentile for AD nuclei
AD.nuclei.genes <- rownames(AD.nuclei.percentile)

# Save the genenames that are in the 25th percentile for AD cells
AD.cells.genes <- rownames(AD.cells.percentile)

# Save the genenames that are in the 25th percentile for WT nuclei
WT.nuclei.genes <- rownames(WT.nuclei.percentile)

# Save the genenames that are in the 25th percentile for WT cells
WT.cells.genes <- rownames(WT.cells.percentile)

# Put all genenames from the percentiles in a list
geneList <- list(AD.nuclei.genes, AD.cells.genes, WT.nuclei.genes, WT.cells.genes)

# Rename the vectors
names(geneList) <- c("AD.nuclei", "AD.cells", "WT.nuclei", "WT.cells")

# Create a venn diagram with the genenames of each group
venn.plot <- venn.diagram(geneList , NULL, cex = 2, category.names=c("AD.nuclei", "AD.cells", "WT.nuclei", "WT.cells"), main="Overlap of the 25th percentile genes in each group", main.cex = 1.5)

# Plot the venn diagram
grid.draw(venn.plot)

# Get the list of genes present in each Venn compartment with the use of the gplot package
genes <- venn(geneList, show.plot=FALSE)

# Save the intersections of the genes object
inters <- attr(genes, "intersections")








