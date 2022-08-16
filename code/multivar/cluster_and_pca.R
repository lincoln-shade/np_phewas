#==============================
# get list of phenotypes with
# sample sizes >5000
# for genetic correlation
#==============================

library(data.table)
library(magrittr)
library(psych)

bin_pheno <- fread("data/mega/mega_np.pheno", na.strings = "-1")
ord_pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = "-1")
ord_pheno[, `:=`(vcid  = NULL, part_def = NULL)] # made from other variables
pheno <- merge(ord_pheno, bin_pheno[, .(FID, IID, hs, grossinf, microinf)],
               by = c('FID', 'IID'))
covar <- fread("data/mega/mega_np.covar")

# PCA on all variables except LATE-NC (has a lot of missingess)
pheno_pca <- pheno[, -c("late", "FID", "IID")]

pca <- principal(pheno_pca, 
                cor = "mixed",
                nfactors = ncol(pheno_pca),
                rotate = "none") 

# examine dendrogram of phenotypic correlations to identify clusters.
# Two major clusters found corresponding roughly to:
#   1. Alzheimer's NP : Braak, CERAD, diffuse amyloid, and CAA
#   2. Vascular NP    : gross/micro infarcts, athero/arteriolosclerosis
pheno_dist <- dist(polychoric(pheno_pca)$rho)
plot(hclust(pheno_dist))

# PCA on Alzheimer's cluster variables
pheno_adnp <- pheno_pca[, .(braak, cerad, caa_ord, diffuse_abeta)]

pca_adnp <- principal(pheno_adnp, 
                 cor = "mixed",
                 nfactors = ncol(pheno_adnp),
                 rotate = "none")

pcs_adnp <- as.data.table(pca_adnp$scores)
pcs_adnp <- cbind(pheno[, .(FID, IID)], pcs_adnp)
pcs_adnp <- merge(pcs_adnp, covar, c("FID", "IID"))
fwrite(pcs_adnp[complete.cases(pcs_adnp)], 
       file = "data/mega/multivar/pcs_adnp.txt", 
       sep = " ", 
       quote = FALSE,
       na = "NA")

# PCA on vascular NP cluster variables
# (leave out arteriolosclerosis because of missingness)
pheno_vcid <- pheno_pca[, .(atheroscler, grossinf, microinf, arteriol)]

pca_vcid <- principal(pheno_vcid, 
                      cor = "mixed",
                      nfactors = ncol(pheno_vcid),
                      rotate = "none") 
pcs_vcid <- as.data.table(pca_vcid$scores)
pcs_vcid <- cbind(pheno[, .(FID, IID)],  pcs_vcid)
pcs_vcid <- merge(pcs_vcid, covar, c("FID", "IID"))
fwrite(pcs_vcid[complete.cases(pcs_vcid)], 
       file = "data/mega/multivar/pcs_vcid.txt", 
       sep = " ", 
       quote = FALSE,
       na = "NA")
