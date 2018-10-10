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
bioPkgTest("scater")

########################################################
#                    Normalization                     #
########################################################

# Load sceset without the quality control assays since they are not neccesary for further analysis
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

# Lognormalize the raw counts and add as assay to the sceset
assay(sceset, "logcounts") <- log2(counts(sceset) + 1)

# Calculate counts per million (CPM) and add as assay to the sceset
assay(sceset, "CPM") <- log2(calculateCPM(sceset) + 1)

# Explore which normalization is best using a PCA plot
plotPCA(sceset, exprs_values = "CPM", colour_by = "Type")

plotPCA(sceset, exprs_values = "logcounts", colour_by = "Type")

# Remove logcounts assay since CPM is best in this case
assay(sceset, "logcounts") <- NULL

# Save sce object for further analysis in another script
saveRDS(sceset, "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_Normalized.rds")
