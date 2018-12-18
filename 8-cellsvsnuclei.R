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
bioPkgTest("edgeR")
bioPkgTest("ComplexHeatmap")


########################################################
#              Celltype Identification                 #
########################################################

# Load quality control non-normalized single cell experiment object
sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE_QC.rds")

# Normalize using RPKM based on the total amount of genes
x_rpkm <- rpkm(assay(sceset, "counts"), 28692)
rpkmsum <- sum(x_rpkm, na.rm=F)
tpm.values <- x_rpkm/rpkmsum * 10^6
# Set RPKM as assay in the sce set
assay(sceset, "RPKM") <- tpm.values

# Remove the RPKM data
rm(x_rpkm)
rm(tpm.values)
gc()

# Create an WT only dataset
WT <- filter(sceset, Condition == "WT")

# Remove the single cell experiment object
rm("sceset")
gc()

# Plot the top 50 genes that are highly expressed
plotHighestExprs(WT)

# Put the top 50 genes in a vector (genes of interest)
goi <- c("Malat1", "Cst3", "Apoe", "Ctsd", "Hexb", "C1qa", "Ctss", "C1qb", "C1qc", "Tmsb4x", "Itm2b", "Lgmn", "Ctsb", "B2m", "Trem2", "Cd81", "Tyrobp", "Ctsz", "Cd9", "Rps29", "Fth1", "Ft1", "mt.Co3", "Ctsl", "mt.Atp6", "Psap", "Foer1g", "Selplg", "H2.D1", "Lyz2", "Rpl41", "Foris", "Cd68", "Csf1r", "Serpine2", "Ly86", "Cst7", "Laptm5", "Grn", "Rps27", "P2ry12", "Rpl18a", "Cd63", "Actb", "Sparc", "Rpl37a", "Hexa", "Cd74", "mt.Co1", "Olfml3")

# Create a logical which genes to use 
genes <- rownames(WT) %in% goi

# Keep the genes of interest in the sceset
WT <- WT[genes, ]

# Filter cells and nuclei in different scesets
cells <- filter(WT, Type == "Microglia")

nuclei <- filter(WT, Type == "Nuclei")

# Create a heatmap of the top50 genes for nuclei and cells to compare
Heatmap(assay(cells, "RPKM"), show_column_names = FALSE) + Heatmap(assay(nuclei, "RPKM"), show_column_names = FALSE)

