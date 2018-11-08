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

########################################################
#                Creating SCE object                   #
########################################################

# Location of the barcodes, genes and matrix files (output from cellranger)
location <- "/media/imgorter/BD_1T/Iris/Data_abbvie/AD_WT/outs/filtered_gene_bc_matrices_mex/mm10"

# Load the cellbarcodes, genenames and molecules files  
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

# Save object to import in another script
saveRDS(sceset, file = "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

# overwrite the rownames to the ensembl molecules (first column)
rownames(molecules) <- make.names(genenames[,1], unique = TRUE)

sceset <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = anno)

# Save as a different object
saveRDS(sceset, file = "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_ensembl.rds")
