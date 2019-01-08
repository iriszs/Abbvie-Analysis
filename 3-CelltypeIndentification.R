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
bioPkgTest("GSEABase")
bioPkgTest("AUCell")
bioPkgTest("doMC")
bioPkgTest("doRNG")

########################################################
#              Celltype Identification                 #
########################################################

sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

# Location of the gene matrix file
gmtFile <- "/media/imgorter/BD_1T/Iris/Scripts/microglia_signatures.gmt"
# Load the gene matrix
geneSets <- getGmt(gmtFile)

# Create a collection of genes that only exist in the dataframe
geneSets <- subsetGeneSets(geneSets, rownames(sceset))
cbind(nGenes(geneSets))

# Build expression based rankings for all genes in each cell
cells_rankings <- AUCell_buildRankings(sceset, nCores=1, plotStats=FALSE)

# Calculate the 'Area under the curve'
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# Calculate thresholds for each gene set based on the AUC
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)

# Explore thresholds  
cells_assignment$Microglia_lavin$aucThr$thresholds

# Cellbarcodes based on the chosen threshold. These are microglia.
selectedCells <- names(which(getAUC(cells_AUC)["Microglia_lavin",]>0.02824006))

# calculate how many cells are not assigned as microglia
length(colnames(sceset)) - length(selectedCells)





