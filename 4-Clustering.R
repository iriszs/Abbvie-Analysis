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
bioPkgTest("Seurat")

########################################################
#                     Clustering                       #
########################################################

# Load normalized single cell experiment (SCE) object
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_Normalized.rds")

# Convert SCE object to Seurat object
seuset <- CreateSeuratObject(assay(sceset, "CPM"), project = "Abbvie")

# Remove SCE object
rm("sceset")
# Collect garbage to clear memory
gc()

# Scale and center data
# Gives warning that the data has not been normalized.
# The data has been normalized however, just not with the seurat function
seuset <- ScaleData(object = seuset)

# Define variable genes
seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR)

length(seuset@var.genes)

# Run PCA dimensionality reduction and save in the Seurat object
seuset <- RunPCA(object = seuset, do.print = FALSE)

# Identify clusters of cells by a shared nearest neighbor modularity optimization based clustering algorithm
# Use the first 4 dimensions
seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:4, resolution = 1.0, save.SNN = TRUE, print.output = 0)

# Print parameters used for the FindClusters function
PrintFindClustersParams(object = seuset)
table(seuset@ident)

# Run TSNE with the first 4 dimensions
seuset <- RunTSNE(object = seuset, dims.use = 1:4, do.fast = TRUE, check_duplicates = FALSE)

# Save Seurat object with the cluster information
saveRDS(seuset, "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/seuset_clustered.rds")

