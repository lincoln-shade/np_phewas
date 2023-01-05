#==================================================
# Figure S1: Heatmap of NPE correlations
#==================================================

library(data.table)
library(magrittr)
library(psych)
library(pheatmap)

bin_pheno <- fread("data/mega/mega_np.pheno", na.strings = "-1")
ord_pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = "-1")
pheno <- merge(ord_pheno, bin_pheno[, .(FID, IID, hs, grossinf, microinf)],
               by = c('FID', 'IID'))
covar <- fread("data/mega/mega_np.covar")

# PCA on all variables except LATE-NC (has a lot of missingess)
pheno_pca <- pheno[, -c("late", "FID", "IID")]
labels = c("Braak NFT Stage", "CERAD Score", "Diffuse Amyloid Plaques", 
           "Arteriolosclerosis",  "Atherosclerosis", "Lewy Body", 
           "CAA", "HS", "Gross Infarction", "Microinfarct")
setnames(pheno_pca, colnames(pheno_pca), labels)

pca <- principal(pheno_pca, 
                 cor = "mixed",
                 nfactors = ncol(pheno_pca),
                 rotate = "none") 

# examine dendrogram of phenotypic correlations to identify clusters.
# Two major clusters found corresponding roughly to:
#   1. Alzheimer's NP : Braak, CERAD, diffuse amyloid, and CAA
#   2. Vascular NP    : gross/micro infarcts, athero/arteriolosclerosis

dendrogram <- pheatmap(polychoric(pheno_pca)$rho, 
                       # labels_row = labels,
                       # labels_col = rep("", 10),
                       filename = "./doc/figure_S1_dendrogram.png",
                       height = 7,
                       width = 8)

dendrogram
