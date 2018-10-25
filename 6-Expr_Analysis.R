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


########################################################
#                Expression Analysis                   #
########################################################

# Load non-normalized single cell experiment object
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

# Remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(sceset) > 0) > 0
sceset <- sceset[keep_feature, ]


# Put AD and WT in different datasets
AD <- filter(sceset, Condition == "AD")
WT <- filter(sceset, Condition == "WT")

# Remove the single cell experiment object
rm("sceset")
gc()

# Create a dataframe for each cell type and condition 
AD.nuclei <- filter(AD, Type == "Nuclei")
AD.cells <- filter(AD, Type == "Microglia")

WT.nuclei <- filter(WT, Type == "Nuclei")
WT.cells <- filter(WT, Type == "Microglia")

# Calculate the mean of the expression of each gene 
AD.nuclei_mean <- rowMeans(assay(AD.nuclei))
# Convert the set to a dataframe
AD.nuclei_mean <- as.data.frame(AD.nuclei_mean)
# Log-normalize the mean
AD.nuclei_mean <- log2(AD.nuclei_mean + 1)
# Set the columnames of the dataframe
colnames(AD.nuclei_mean) <- "AD.nuclei"

# Calculate the mean of the expression of each gene 
AD.cells_mean <- rowMeans(assay(AD.cells))
# Convert the set to a dataframe
AD.cells_mean <- as.data.frame(AD.nuclei_mean)
# Log-normalize the mean
AD.cells_mean <- log2(AD.cells_mean + 1)
# Set the columnames of the dataframe
colnames(AD.cells_mean) <- "AD.cells"

# Calculate the mean of the expression of each gene
WT.nuclei_mean <- rowMeans(assay(WT.nuclei))
# Convert the set to a dataframe
WT.nuclei_mean <- as.data.frame(WT.nuclei_mean)
# Log-normalize the mean
WT.nuclei_mean <- log2(WT.nuclei_mean + 1)
# Set the columnames of the dataframe
colnames(WT.nuclei_mean) <- "WT.nuclei"

# Calculate the mean of the expression of each gene
WT.cells_mean <- rowMeans(assay(WT.cells))
# Convert the set to a dataframe
WT.cells_mean <- as.data.frame(WT.cells_mean)
# Log-normalize the mean
WT.cells_mean <- log2(WT.cells_mean + 1)
# Set the columnames of the dataframe
colnames(WT.cells_mean) <- "WT.cells"

# Combine the normalized mean for each condition and celltype together
df <- cbind(AD.nuclei_mean, AD.cells_mean, WT.nuclei_mean, WT.cells_mean)

# Check if all columns are numeric
sapply(df, class)

# Plot the nuclei vs cells in AD
ggplot(df, aes(x = df$AD.nuclei, y = df$AD.cells)) + geom_point(color = "blue", shape = 1) + geom_smooth(method = "lm", color = "red") + ggtitle("Log-mean expression of nuclei and cells in AD") + ylab("AD cells") + xlab("AD nuclei")

# Plot the nuclei vs cells in WT
ggplot(df, aes(x = df$WT.nuclei, y = df$WT.cells)) + geom_point(color = "blue", shape = 1) + geom_smooth(method = "lm", color = "red") + ggtitle("Log-mean expression of nuclei and cells in WT") + ylab("WT cells") + xlab("WT nuclei")
