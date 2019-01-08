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
bioPkgTest("ggplot2")
bioPkgTest("plotly")


########################################################
#                Expression Analysis                   #
########################################################

# Load non-normalized single cell experiment object
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

# Remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(sceset) > 0) > 0
sceset <- sceset[keep_feature, ]

# Put APP and WT in different datasets
APP <- scater::filter(sceset, Condition == "APP")
WT <- scater::filter(sceset, Condition == "WT")

# Remove the single cell experiment object
rm("sceset")
gc()

# Create a sceset for each cell type and condition 
APP.nuclei <- scater::filter(APP, Type == "Nuclei")
APP.cells <- scater::filter(APP, Type == "Cells")

WT.nuclei <- scater::filter(WT, Type == "Nuclei")
WT.cells <- scater::filter(WT, Type == "Cells")

# Put all scesets in a list
types <- list(APP.nuclei, APP.cells, WT.nuclei, WT.cells)

# Create a matrix with the size of the four scesets and the number of genes that is left
df <- data.frame(matrix(ncol = 4, nrow = 14245))

# For each sceset
for (i in 1:4){
  # Calculate the rowmeans of the counts
  mean <- rowMeans(assay(types[[i]], "counts"))
  # Set as dataframe
  mean <- as.data.frame(mean)
  # Log-normalize
  mean <- log2(mean +1)
  # Add to the previously made dataframe
  df[i] <- mean
}

# Set colnames of dataframe
colnames(df) <- c("APP.nuclei", "APP.cells", "WT.nuclei", "WT.cells")

# Set the rownames of the last made mean as gene in the dataframe to make plot interactive
df$gene <- rownames(mean)

# Plot the nuclei vs cells in AD
ggplotly(ggplot(df, aes(x = df$APP.nuclei, y = df$APP.cells, label = df$gene)) + geom_point(color = "blue", shape = 1) + geom_abline(color = "red") + ggtitle("Log-mean expression of nuclei and cells in APP") + ylab("APP cells") + xlab("APP nuclei") + xlim(0, 3) + ylim(0, 3))

# Plot the nuclei vs cells in WT
ggplotly(ggplot(df, aes(x = df$WT.nuclei, y = df$WT.cells, label = df$gene)) + geom_point(color = "blue", shape = 1) + geom_abline(color = "red") + ggtitle("Log-mean expression of nuclei and cells in WT") + ylab("WT cells") + xlab("WT nuclei") + xlim(0, 6) + ylim(0, 6))

