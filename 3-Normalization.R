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
bioPkgTest("edgeR")
bioPkgTest("scater")

########################################################
#                    Normalization                     #
########################################################

# Load sceset without the quality control assays since they are not neccesary for further analysis
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

x_rpkm <- rpkm(assay(sceset, "counts"), 28692)
rpkmsum <- sum(x_rpkm, na.rm=F)
tpm.values <- x_rpkm/rpkmsum * 10^6
#tpm.values <- tpm.values[,order(colnames(tpm.values))]

assay(sceset, "RPKM") <- tpm.values

plotPCA(sceset, exprs_values = "RPKM", colour_by = "Type")

# Save sce object for further analysis in another script
saveRDS(sceset, "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_Normalized.rds")


sceset <- readRDS
