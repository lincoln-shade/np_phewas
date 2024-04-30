#==================================================
# Figure S1: Heatmap of NPE correlations
#==================================================

library(data.table)
library(magrittr)
library(psych)
library(pheatmap)
library(polycor)

pheno <- fread("shared/nacc_rosmap_act_np.pheno", na.strings = "-1")
covar <- fread("data/mega/mega_np.covar")

# PCA on all variables
pheno_pca <- pheno[, -c("FID", "IID")]
labels = c("LATE-NC", "Braak NFT Stage", "CERAD Score", "Amyloid-Beta Plaques", 
           "Arteriolosclerosis","Atherosclerosis", "Lewy Body", 
           "CAA", "HS", "Microinfarct", "Gross Infarction")
setnames(pheno_pca, colnames(pheno_pca), labels)

pca <- principal(pheno_pca, 
                 cor = "mixed",
                 nfactors = ncol(pheno_pca),
                 rotate = "none") 

# examine dendrogram of phenotypic correlations to identify clusters.
# Two major clusters found corresponding roughly to:
#   1. Alzheimer's NP : Braak, CERAD, diffuse amyloid, and CAA
#   2. Vascular NP    : gross/micro infarcts, athero/arteriolosclerosis
cors <- polychoric(pheno_pca)
dendrogram <- pheatmap(cors$rho, 
                       # labels_row = labels,
                       # labels_col = rep("", 10),
                       filename = "./doc/figure_S1_dendrogram.png",
                       height = 7,
                       width = 8, 
                       angle_col = "45")

dendrogram

# test specific pairs for significance
test_cor <- function(v1, v2) {
  cors <- polychor(v1, v2, ML=TRUE, std.err = TRUE)
  pval <- pchisq(cors$chisq, cors$df, lower.tail = F)
  return(c(cors$rho, pval))
}
braak_microinf <- test_cor(pheno_pca$`Braak NFT Stage`, pheno_pca$Microinfarct)

# compare correlations across data sources
pheno_nacc <- pheno[FID == "0"]
pheno_rosmap <- pheno[FID == "1"]
pheno_act <- pheno[FID == "2"]

late_hs <- test_cor(pheno$late, pheno$hs)
late_hs_nacc <- test_cor(pheno_nacc$late, pheno_nacc$hs)
late_hs_rosmap <- test_cor(pheno_rosmap$late, pheno_rosmap$hs)
late_hs_act = test_cor(pheno_act$late, pheno_act$hs)

pheno[, test_cor(braak, atheroscler)]
