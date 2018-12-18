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

cells <- c("APP1-cells", "APP2-cells")
nuclei <- c("APP1-nuclei", "APP2-nuclei")

DAM_analyse <- function(samples, type){
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
  
  # Determine mitochondrial genes
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuset@data), value = TRUE)
  # Calculate percentage of mitochondrial genes
  percent.mito <- Matrix::colSums(seuset@raw.data[mito.genes, ])/Matrix::colSums(seuset@raw.data)
  
  # Add mitochondrial percentage to meta data
  seuset <- AddMetaData(object = seuset, metadata = percent.mito, col.name = "percent.mito")
  
  # multiply the expression by 10.000 and lognormalize.
  seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find highly variable genes
  # No cutoffs used
  # Don't add to Seurat object, we want to use all genes for the next steps.
  var_genes <- FindVariableGenes(object = seuset, mean.function = ExpMean, dispersion.function = LogVMR, set.var.genes = FALSE)
  
  # Regress the data on number of detected molecules per cell and the mitochondrial gene content percentage
  seuset <- ScaleData(object = seuset, vars.to.regress = c("nUMI", "percent.mito"))
  
  #Sys.sleep(180)
  
  ########################################################
  #        Perform linear dimensional reduction          #
  ########################################################
  
  pdf(paste0("/media/imgorter/BD_1T/Iris/Results/DAM_analysis/workflow//", type, ".pdf"))
  
  # Run Princincipal component analysis on all genes
  seuset <- RunPCA(object = seuset, pc.genes = rownames(seuset@data), do.print = FALSE)
  
  # Visualize the genes associated with first two principal components
  VizPCA(object = seuset, pcs.use = 1:2)
  
  # Create PCA plot with the first two principal components
  PCAPlot(object = seuset, dim.1 = 1, dim.2 = 2)
  
  # runs the pre-computed PCA to the entire dataset
  seuset <- ProjectPCA(object = seuset, do.print = FALSE)
  
  # Create heatmap with the first two princripal components and use the first 500 cells
  PCHeatmap(object = seuset, pc.use = 1:2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
  
  # Create a PCElbow plot to identify which dimensions matter
  PCElbowPlot(object = seuset)
  
  ########################################################
  #                    Clustering                        #
  ########################################################
  
  # Identify cluster of cells by a SNN modularity optimization
  # Use the first 15 dimensions (determined from elbow plot)
  # Low resolution to get a smaller number of communities
  seuset <- FindClusters(object = seuset, reduction.type = "pca", dims.use = 1:15, resolution = 0.1, print.output = 0, save.SNN = TRUE)
  
  # Run tsne dimensionality reduction on all genes
  # Use the first 15 dimensions
  seuset <- RunTSNE(object = seuset, dims.use = 1:15, do.fast = TRUE, check_duplicates = FALSE)
  
  # Cluster based tsne plot
  TSNEPlot(object = seuset)
  
  # Genes that are suppressed in neurodegenerative diseases (Butovsky et al.)
  suppressed <- c("P2ry12", "Ccr5", "Cd33", "Csf1r", "Cx3cr1", "Glul",  "Gpr34", "Adgrg1", "Tgfb1", "Tgfbr1", "Serinc3", "Siglech", "Mertk", "Bin1", "Tmem119", "Pu.1", "Sall1", "Mafb", "Smad3", "Mef2a", "Egr1", "Jun")
  
  # Create heatmap of the suppressed genes
  hm1 <- DoHeatmap(seuset, genes.use = suppressed, slim.col.label = TRUE, group.spacing = 0.5, col.low = "Red", col.mid = "Black", col.high = "Green", title = paste0("Heatmap of APP microglia ", type, " with suppressed genes"))
  print(hm1)
  
  # Genes that are induced in neurodegenerative diseases (Butovsky et al.)
  induced <- c("Apoe", "Axl", "Bhlhe40", "Clec7a", "Csf1", "Cst7", "Ctsb", "Ctsd", "Ctsl", "Cybb", "Fabp5", "Fth1", "Itgax", "Gnas", "Gpnmb", "Grn", "Il1b", "Lgals3", "Lilrb4", "Lpl", "Lyz2", "Mir155", "Msr1", "Nos2", "Spp1", "Tfec", "Trem2", "Tyrobp", "Vegfa")
  
  # Create heatmap of induced genes
  hm2 <- DoHeatmap(seuset, genes.use = induced, slim.col.label = TRUE, group.spacing = 0.5, col.low = "Red", col.mid = "Black", col.high = "Green", title = paste0("Heatmap of APP microglia ", type,  " with induced genes"))
  print(hm2)
  
  # Homeostatic microglia genes (Keren-Shaul et al.)
  homeostatic <- c("Hexb", "Cst3", "Cx3cr1", "Ctsd", "Csf1r", "Ctss", "Sparc", "Tmsb4x", "P2ry12", "C1qa", "C1qb")
  
  # Create heatmap of homeostatic microglia genes
  hm3 <- DoHeatmap(seuset, genes.use = homeostatic, slim.col.label = TRUE, group.spacing = 0.5, col.low = "Red", col.mid = "Black", col.high = "Green", title = paste0("Heatmap of APP microglia ", type, " of homeostatic mg genes"))
  print(hm3)
  
  # Genes that are involved (higher expression) in stage 2 Disease Associated Microglia (DAM) (Keren-Shaul)
  stage2 <- c("Trem2", "Axl", "Cts7", "Ctsl", "Lpl", "Cd9", "Csf1", "Ccl6", "Itgax", "Clec7", "Lilrb4", "Timp2")
  
  # Create heatmap of the stage 2 DAM genes
  hm4 <- DoHeatmap(seuset, genes.use = stage2, slim.col.label = TRUE, group.spacing = 0.5, col.low = "Red", col.mid = "Black", col.high = "Green", title = paste0("Heatmap of APP microglia ", type, " with stage 2 DAM genes"))
  print(hm4)
  
  # Find variable genes between clusters
  markers <- FindAllMarkers(seuset)
  
  saveRDS(seuset, file = paste0("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/", type, "_DAM.rds" ))
  
  dev.off()
  
}

# Execute the analysis for both cells and nuclei
DAM_analyse(cells, "cells")
DAM_analyse(nuclei, "nuclei")



