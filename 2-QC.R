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
#                  Quality Control                     #
########################################################

# Load the SCE object from the first script
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

#  Define mitochondrial genes
isSpike(sceset, "mt") <- grep(pattern = "^mt\\.", rownames(sceset), value = TRUE)

# Calculate QC for mouse mitochondrial genes
sceset <- calculateQCMetrics(sceset, feature_controls = list(MT = isSpike(sceset, "mt")))

# Calculate quality metrics for each cell and feature
sceset <- calculateQCMetrics(sceset)

# Create histogram of the RNA molecules distribution across cells
hist(sceset$total_counts, breaks = 100, main = "Library Size", xlab = "Counts")

# Create histogram of the unique genes distribution across cells
hist(sceset$total_features_by_counts , breaks = 100, main = "Unique genes detected per cell", xlab = "Number of genes")

median(sceset$total_features_by_counts)

plotColData(sceset, x = "total_features_by_counts", y = "pct_counts_MT", colour_by = "Type") + ggtitle("Percentage of counts in MT genes")

# Filters based on RNA molecules and unique genes
# are 'off' for now, but not removed to make it easy to change it fast
filter_by_total_counts <- (sceset$total_counts > 0)

filter_by_expr_features <- (sceset$total_features > 0)

# Only use cells that pass the filter
sceset$use <- (filter_by_expr_features & filter_by_total_counts)

# Save object to import in another script
saveRDS(sceset, "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_QC.rds" )
