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
bioPkgTest("GeneOverlap")
bioPkgTest("ggplot2")


########################################################
#                        DEA                           #
########################################################

sceset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/SCE.rds")

assay(sceset)[1:5,1:5]

# Put AD and WT in different datasets
AD <- filter(sceset, Condition == "AD")
WT <- filter(sceset, Condition == "WT")

rm("sceset")
gc()

AD.nuclei <- filter(AD, Type == "Nuclei")
AD.cells <- filter(AD, Type == "Microglia")

WT.nuclei <- filter(WT, Type == "Nuclei")
WT.cells <- filter(WT, Type == "Microglia")

df <- cbind(assay(AD.nuclei), assay(AD.cells), assay(WT.nuclei), assay(WT.cells))

apply(df, MARGIN = c(1, 2), FUN = mean)

keep_feature <- rowSums(counts(AD.nuclei) > 0) > 0
AD.nuclei <- AD.nuclei[keep_feature, ]

keep_feature <- rowSums(counts(AD.cells) > 0) > 0
AD.cells <- AD.cells[keep_feature, ]

keep_feature <- rowSums(counts(WT.nuclei) > 0) > 0
WT.nuclei <- WT.nuclei[keep_feature, ]

keep_feature <- rowSums(counts(WT.cells) > 0) > 0
WT.cells <- WT.cells[keep_feature, ]

overlapAD.obj <- newGeneOverlap(rownames(AD.nuclei), rownames(AD.cells), genome.size = 28692)

overlapAD.obj <- testGeneOverlap(overlapAD.obj)

overlapAD <- overlapAD.obj@intersection

overlapWT.obj <- newGeneOverlap(rownames(WT.nuclei), rownames(WT.cells), genome.size = 28692)

overlapWT.obj <- testGeneOverlap(overlapWT.obj)

overlapWT <- overlapWT.obj@intersection

WT.cells.NO <- length(rownames(WT.cells)) - length(overlapWT)

WT.nuclei.NO <- length(rownames(WT.nuclei)) - length(overlapWT)

AD.cells.NO <- length(rownames(AD.cells)) - length(overlapAD)

AD.nuclei.NO <- length(rownames(AD.nuclei)) - length(overlapAD)

length(overlapAD)

as.data.frame(counts(filter(AD, Type == "Microglia")))

AD.plot <- ggplot(as.data.frame(counts(AD)), aes(x = as.data.frame(counts(AD.cells)), y = as.data.frame(counts(AD.nuclei)))) + geom_point()

AD.plot


