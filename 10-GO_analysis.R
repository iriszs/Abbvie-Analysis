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
bioPkgTest("Seurat")
bioPkgTest("clusterProfiler")
bioPkgTest("org.Mm.eg.db")
bioPkgTest("edgeR")
bioPkgTest("scater")
pkgTest("xlsx")

########################################################
#                      Functions                       #
########################################################

GO <- function(counts){
  # Calculate the 25th percentile 
  percentile <- counts[rowSums(counts >0) > quantile(rowSums(counts) , probs=0.75),]
  
  # Get the rownames (which are the gene symbols) of the 25th percentile
  symbol <- rownames(percentile)
  
  # Put those genes in a dataframe
  genes.df <- as.data.frame(symbol)
  
  # Convert the gene symbols to entrez id's and add to the genes dataframe
  genes.df[,"entrez"] <- mapIds(org.Mm.eg.db, keys = genes.df$symbol, column="ENTREZID", keytype = "SYMBOL", multiVals="first")
  
  # Molecular function
  GOEA <- enrichGO(as.character(genes.df$entrez), OrgDb="org.Mm.eg.db", ont = "MF")
  
  # Plot the results in a dotplot
  print(dotplot(GOEA, title = "GOEA of Molecular Function"))
  
  # Biological Process
  GOEA <- enrichGO(as.character(genes.df$entrez), OrgDb="org.Mm.eg.db", ont = "BP")
  
  # Plot the results in a dotplot
  print(dotplot(GOEA, title = "GOEA of Biological Processes"))
  
  # Cellular Component
  GOEA <- enrichGO(as.character(genes.df$entrez), OrgDb="org.Mm.eg.db", ont = "CC")
  
  # Plot the results in a dotplot
  print(dotplot(GOEA, title = "GOEA of Cellular Component"))
}

########################################################
#                     GO Analysis                      #
########################################################


########################################################
#                      All Genes                       #
########################################################
# Load RDS
seuset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/cells_DAM.rds")

# Get all cells that are assigned to cluster 0
damSet <- SubsetData(seuset, ident.use = 0)

# Calculate RPKM values for the counts based on the total number of genes
x_rpkm <- rpkm(damSet@raw.data, 28692)
rpkmsum <- sum(x_rpkm, na.rm=F)
tpm.values <- x_rpkm/rpkmsum * 10^6

# Set the RPKM values as the data assay
damSet@data <- tpm.values

# Call GO function as defined above
GO(damSet@data)

# Calculate the mean expression of each gene
mean <- rowMeans(damSet@data)
mean <- as.data.frame(mean)
# Order the expression from highest to lowest
mean <- mean[order(-mean$mean), , drop = FALSE]

# Write the 3000 genes that have the average highest expression to an xlsx file to use in Metascape
write.xlsx(rownames(mean)[0:3000], file = "/media/imgorter/BD_1T/Iris/Results/DAM_analysis/GO_analysis/first_3000.xlsx", row.names = FALSE, col.names = FALSE)

########################################################
#            Without ribo and mito genes               #
########################################################

# Convert the RPKM values and genes to a single cell experiment object
# Since it is not possible to remove genes from a Seurat object
sceset <- SingleCellExperiment(assays = list(counts = damSet@data))

# Find the mitochondrial genes (those genes start with mt-)
is.mito <- grepl("^mt-", rownames(sceset))

# Remove those genes from the sceset
sceset <- sceset[is.mito==FALSE, ]

# Find the ribosomal genes(those genes start with Rp and have a number in them)
is.ribo <- grepl("^Rp[sl][[:digit:]]", rownames(sceset))

# Remove the ribosomal genes from the sceset
sceset <- sceset[is.ribo==FALSE, ]

# Calculate how many genes are removed
length(which(is.mito)) + length(which(is.ribo))

# Call GO function as defined above
GO(assay(sceset, "counts"))

# Calculate the mean expression of each gene
mean <- rowMeans(assay(sceset, "counts"))
mean <- as.data.frame(mean)
# Order the expression from highest to lowest
mean <- mean[order(-mean$mean), , drop = FALSE]

# Write the 3000 genes that have the average highest expression to an xlsx file to use in Metascape
write.xlsx(rownames(mean)[0:3000], file = "/media/imgorter/BD_1T/Iris/Results/DAM_analysis/GO_analysis/first_3000_without_ribo_mito.xlsx", row.names = FALSE, col.names = FALSE)

