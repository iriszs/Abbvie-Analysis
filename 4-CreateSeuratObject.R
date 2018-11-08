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

########################################################
#                 Setup Seurat Object                  #
########################################################

# location of dataset
dataset_loc <- "/media/imgorter/BD_1T/Iris/Data_abbvie/mouse_sc_microglia_nuclei"
# Non-aggregated sample names
samples <- c("WT1-cells", "WT1-nuclei", "WT2-cells", "WT2-nuclei", "APP1-cells", "APP1-nuclei", "APP2-cells", "APP2-nuclei")

# Source of loading in data: https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART1.html
data <- sapply(samples, function(i){
  d10x <- Read10X(file.path(dataset_loc, i, "outs/filtered_gene_bc_matrices/mm10/"))
  # Add samplename to the colnames
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

# Combine the datasets together
seuset.data <- do.call("cbind", data)

# Create Seurat object
seuset <- CreateSeuratObject(
  seuset.data,
  project = "Abbvie",
  names.field = 2,
  names.delim = "\\-")

# Read in metadata
meta.data <- read.csv("/media/imgorter/BD_1T/Iris/Data_abbvie/reference_AD_WT.csv", header = TRUE)

# For each column in the metadata file
for (j in 1:ncol(meta.data)){
  # Create a dataframe with the lenght of the data in the Seurat object
  meta <- data.frame(matrix(NA, nrow = length(colnames(seuset@raw.data)), ncol = 1))
  # Add all meta data rows to the dataframe
  t <- meta.data[names(samples[1]),]
  variable <- as.character(t[,j])
  vector <- c(rep(variable, length(seuset@cell.names)))
  meta[,1] <- as.factor(vector)
  colnames(meta) <- c(as.character(colnames(t[j])))
  rownames(meta) <- colnames(seuset@raw.data)
  seuset <- AddMetaData(seuset, metadata = meta)
}

# Save Seurat object
saveRDS(seuset, file = "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/Seurat.rds")
