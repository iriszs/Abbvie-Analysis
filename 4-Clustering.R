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
bioPkgTest("")

########################################################
#                     Clustering                       #
########################################################

sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_Normalized.rds")

seuset <- CreateSeuratObject(assay(sceset, "CPM"), project = "Abbvie")
rm("sceset")
gc()

# Gives warning that the data has not been normalized.
# The data has been normalized however, just not with the seurat function
seuset <- ScaleData(object = seuset)
#regress out proberen om effect van aantal genen per cel er uit te filteren. 

seuset <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR)

FindVariableGenes(seuset, set.var.genes = FALSE)

length(seuset@var.genes)

seuset <- RunPCA(object = seuset, do.print = FALSE)

seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:4, resolution = 1.0, save.SNN = TRUE, print.output = 0)
PrintFindClustersParams(object = seuset)
table(seuset@ident)

seuset <- RunTSNE(object = seuset, dims.use = 1:4, do.fast = TRUE, check_duplicates = FALSE)
#group by ngene. 
#kleuren op sample

saveRDS(seuset, "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/seuset_clustered.rds")

