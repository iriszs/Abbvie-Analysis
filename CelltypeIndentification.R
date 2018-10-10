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

# Location of data
location <- "/media/imgorter/BD_1T/Iris/Data_abbvie/AD_WT/outs/filtered_gene_bc_matrices_mex/mm10"

# Load the data  
cellbarcodes <- read.table(paste0(location, "/barcodes.tsv"))
genenames <- read.table(paste0(location, "/genes.tsv"))
molecules <- Matrix::readMM(paste0(location, "/matrix.mtx"))

# Load reference/annotation file
anno <- read.csv("/media/imgorter/BD_1T/Iris/Data_abbvie/reference_AD_WT.csv")

# Set the rownames of the molecules matrix as the second column of the genenames matrix (which are the gene symbols)
rownames(molecules) <- make.names(genenames[,2], unique = TRUE)
# Set colnames of molecules as the first column of the cellbarcodes matrix (which are the cellbarcodes)
colnames(molecules) <- cellbarcodes[,1]

# Create sceset object using the molecules matrix and the annotation info
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = anno)

# Location of the gene matrix file
gmtFile <- "/media/imgorter/BD_1T/Iris/microglia_signatures.gmt"
#gmtFile <- paste(file.path(system.file('examples', package='AUCell')), "geneSignatures.gmt", sep="/")
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


# Select cells based on the chosen threshold. These are microlgia
selectedCells <- names(which(getAUC(cells_AUC)["Microglia_lavin",]>0.02824006))

# calculate how many cells are not assigned as microglia
length(colnames(sceset)) - length(selectedCells)





