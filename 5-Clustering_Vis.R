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
bioPkgTest("openxlsx")

########################################################
#                     Visualizing                      #
########################################################

# Load seurat object 
seuset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/seuset_clustered.rds")

# Plot tSNE that is saved in the object
TSNEPlot(object = seuset)

# Open PDF
pdf("/media/imgorter/BD_1T/Iris/Results/Genes_tsne/genes.pdf")

# Plot genes that express in cells but not in nuclei
FeaturePlot(object = seuset, features.plot = c("Tmsb4x", "Irak2", "Rplp2", "Bag3"), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
# Plot mitochondrial genes in four different plots to make it more visible
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt\\.", rownames(seuset@data), value = TRUE)[1:4]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt\\.", rownames(seuset@data), value = TRUE)[5:8]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt\\.", rownames(seuset@data), value = TRUE)[9:12]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
FeaturePlot(object = seuset, features.plot = c(grep(pattern = "^mt\\.", rownames(seuset@data), value = TRUE)[13]), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
# Plot AD genes
FeaturePlot(object = seuset, features.plot = c("Apoe"), cols.use = c("grey", "red"), reduction.use = "tsne", no.legend = FALSE)
# Plot homeostatic microglia genes
FeaturePlot(object = seuset, features.plot = c("Trem2", "P2ry12"), cols.use = c("grey", "red", "green"), reduction.use = "tsne", overlay = TRUE, no.legend=FALSE)

# Close PDF
dev.off()

cluster1.markers <- FindMarkers(object = seuset, ident.1 = 0)

cluster2.markers <-FindMarkers(object = seuset, ident.1 = 1)

cluster3.markers <- FindMarkers(object = seuset, ident.1 = 2)

cluster4.markers <- FindMarkers(object= seuset, ident.1 = 3)

write.xlsx(cluster1.markers, "/media/imgorter/BD_1T/Iris/Results/tsne_markers/cluster1.xlsx", col.names=TRUE, row.names=TRUE)
write.xlsx(cluster2.markers, "/media/imgorter/BD_1T/Iris/Results/tsne_markers/cluster2.xlsx", col.names=TRUE, row.names=TRUE)
write.xlsx(cluster3.markers, "/media/imgorter/BD_1T/Iris/Results/tsne_markers/cluster3.xlsx", col.names=TRUE, row.names=TRUE)
write.xlsx(cluster4.markers, "/media/imgorter/BD_1T/Iris/Results/tsne_markers/cluster4.xlsx", col.names=TRUE, row.names=TRUE)

