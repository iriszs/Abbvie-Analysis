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
bioPkgTest("clusterProfiler")
bioPkgTest("org.Mm.eg.db")
bioPkgTest("edgeR")
bioPkgTest("scater")
pkgTest("xlsx")


########################################################
#                     GO Analysis                      #
########################################################

seuset <- readRDS("/media/imgorter/BD_1T/Iris/Scripts/abbvie/RDS/cells_DAM.rds")

damSet <- SubsetData(seuset, ident.use = 0)

x_rpkm <- rpkm(damSet@raw.data, 28692)
rpkmsum <- sum(x_rpkm, na.rm=F)
tpm.values <- x_rpkm/rpkmsum * 10^6

damSet@data <- tpm.values

percentile <- damSet@data[rowSums(damSet@data >0) > quantile(rowSums(damSet@data) , probs=0.75),]

symbol <- rownames(percentile)

genes.df <- as.data.frame(symbol)

genes.df[,"entrez"] <- mapIds(org.Mm.eg.db, keys = genes.df$symbol, column="ENTREZID", keytype = "SYMBOL", multiVals="first")

GOEA <- enrichGO(as.character(genes.df$entrez), OrgDb="org.Mm.eg.db")

dotplot(GOEA)

