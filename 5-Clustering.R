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
pkgTest("ggplot2")

########################################################
#                   Quality Control                    #
########################################################

# Load Seurat object
seuset <- readRDS(file = "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/Seurat.rds")

# Determine mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuset@data), value = TRUE)
# Calculate percentage of mitochondrial genes
percent.mito <- Matrix::colSums(seuset@raw.data[mito.genes, ])/Matrix::colSums(seuset@raw.data)

# Add mitochondrial percentage to meta data
seuset <- AddMetaData(object = seuset, metadata = percent.mito, col.name = "percent.mito")
# nGene and nUMI are pre-calculated when creating the seuset
png(filename="/media/imgorter/BD_1T/Iris/Results/QualityControl/Seurat/violinQC.png")
VlnPlot(object = seuset, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# Plot quality control plots
par(mfrow = c(1, 2))
GenePlot(object = seuset, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seuset, gene1 = "nUMI", gene2 = "nGene")

# multiply the expression by 10.000 and lognormalize.
seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", scale.factor = 10000)

# Find highly variable genes
# No cutoffs used
# Don't add to Seurat object, we want to use all genes for the next steps.
var_genes <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR, set.var.genes = FALSE)

# Regress the data on number of detected molecules per cell and the mitochondrial gene content percentage
seuset <- ScaleData(object = seuset, vars.to.regress = c("nUMI", "percent.mito"))

########################################################
#        Perform linear dimensional reduction          #
########################################################

# Run Princincipal component analysis on all genes
seuset <- RunPCA(object = seuset, pc.genes = rownames(seuset@data), do.print = FALSE)

# Visualize the genes associated with first two principal components
VizPCA(object = seuset, pcs.use = 1:2)

# Create PCA plot with the first two principal components
png(filename="/media/imgorter/BD_1T/Iris/Results/tsne/PCA.png")
PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
dev.off()

# runs the pre-computed PCA to the entire dataset
seuset <- ProjectPCA(object = seuset, do.print = FALSE)

# Create heatmap with the first two princripal components and use the first 500 cells
PCHeatmap(object = seuset, pc.use = 1:2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = seuset, pc.use = 5:6, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

# Create a PCElbow plot to identify which dimensions matter
png(filename="/media/imgorter/BD_1T/Iris/Results/tsne/PCelbow.png")
PCElbowPlot(object = seuset)
dev.off()

########################################################
#                    Clustering                        #
########################################################

# Identify cluster of cells by a SNN modularity optimization
# Use the first 10 dimensions (determined from elbow plot)
# Low resolution to get a smaller number of communities
seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:10, resolution = 0.1, print.output = 0, save.SNN = TRUE)

# Run tsne dimensionality reduction on all genes
# Use the first 10 dimensions
seuset <- RunTSNE(object = seuset, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)

seuset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/Seurat_clustered.rds")

# Plot tsne plot by mouse
png(filename="/media/imgorter/BD_1T/Iris/Results/tsne/tsne_mouse.png")
TSNEPlot(object = seuset, group.by = "orig.ident")
dev.off()

# Cluster based tsne plot
png(filename="/media/imgorter/BD_1T/Iris/Results/tsne/tsne_.png")
TSNEPlot(object = seuset)
dev.off()

cluster.info <- as.data.frame(table(seuset@ident, seuset@meta.data$orig.ident))

colnames(cluster.info) <- c("Cluster", "Mouse", "Cells")

# Barplot that indicates which mice is assigned to which cluster
ggplot(cluster.info, aes(x = Mouse, y = Cells, fill = Cluster)) + geom_bar(stat="identity")

# Calculate percentage of each cluster out of the total number of cells for each mouse
AD1 <- cluster.info[grep("AD1", cluster.info$Mouse), ]
AD1[, "Percentage"] <- AD1[,3] / sum(AD1[,3]) * 100

AD2 <- cluster.info[grep("AD2", cluster.info$Mouse), ]
AD2[, "Percentage"] <- AD2[,3] / sum(AD2[,3]) * 100

WT1 <- cluster.info[grep("WT1", cluster.info$Mouse), ]
WT1[, "Percentage"] <- WT1[,3] / sum(WT1[,3]) * 100

WT2 <- cluster.info[grep("WT2", cluster.info$Mouse), ]
WT2[, "Percentage"] <- WT2[,3] / sum(WT2[,3]) * 100

# Plot relative barplot
cluster.info <- rbind(AD1, AD2, WT1, WT2)

ggplot(cluster.info, aes(x = Mouse, y = Percentage, fill = Cluster)) + geom_bar(stat = "identity")


pdf("/media/imgorter/BD_1T/Iris/Results/tsne/genes.pdf")

FeaturePlot(object = seuset, features.plot = c("Tmsb4x", "Irak2", "Rplp2", "Bag3"), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
# Plot mitochondrial genes 
# four different plots to make it more visible
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt-.", rownames(seuset@data), value = TRUE)[1:4]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt-.", rownames(seuset@data), value = TRUE)[5:8]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt-.", rownames(seuset@data), value = TRUE)[9:12]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt-.", rownames(seuset@data), value = TRUE)[13]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
# Plot AD genes
FeaturePlot(object = seuset, features.plot = c("Apoe"), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
# Plot homeostatic microglia genes
FeaturePlot(object = seuset, features.plot = c("Trem2", "P2ry12"), cols.use = c("grey", "red", "green"), reduction.use = "tsne", overlay = TRUE, no.legend=FALSE)

dev.off()

# Save Seurat object with the clustering information
saveRDS(seuset, file = "/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/Seurat_clustered.rds")
