#==============================
# get list of phenotypes with
# sample sizes >5000
# for genetic correlation
#==============================

library(data.table)
library(magrittr)
library(polycor)

bin_pheno <- fread("data/mega/mega_np.pheno", na.strings = "-1")
ord_pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = "-1")
ord_pheno[, `:=`(vcid  = NULL, part_def = NULL)]
pheno <- merge(ord_pheno, bin_pheno[, .(FID, IID, hs, grossinf, microinf)],
               by = c('FID', 'IID'))
# pheno[, hs := NULL]

for (i in colnames(pheno[, 3:ncol(pheno)])) {
  if (is.numeric(pheno[[i]])) {
    if (pheno[!(is.na(get(i))), 
                  sum(get(..i))] %in% c(0, .N)) {
      pheno[, (i) := NULL]
    }
  }
}


pheno_4k <- character()
for (i in colnames(pheno)) {
  if (pheno[!(is.na(get(i))), .N] >= 2000) {
    pheno_4k <- c(pheno_4k, i)
  }
}
pheno_4k <- pheno[, ..pheno_4k]
pheno_4k <- pheno_4k[complete.cases(pheno_4k)]

n_cols <- ncol(pheno_4k) - 2L
cor_mat <- matrix(nrow = n_cols, ncol = n_cols)
colnames(cor_mat) <- 
  rownames(cor_mat) <- 
  colnames(pheno_4k[, -c("FID", "IID")])
for (i in 1L:n_cols) {
  pheno1 <- colnames(pheno_4k)[i + 2L]
  for (j in i:n_cols) {
    pheno2 <- colnames(pheno_4k)[j + 2L]
    cor_mat[i, j] <- cor_mat[j, i] <- 
      pheno_4k[, polychor(get(..pheno1), 
                          get(..pheno2))]
  }
}

cor_mat[is.na(cor_mat)] <- 0
cor_mat <- round(cor_mat, 2)

pca1 <- prcomp(pheno_4k[, -c("FID", "IID")], scale. = T, center = T)
pca2 <- princomp(cor_mat)
pcs <- data.table(FID = pheno_4k$FID, 
                  IID = pheno_4k$IID, 
                  pca1$x)

pheno_4k[IID %in% pcs[PC1 > 4, IID]]

pheno_dist <- dist(cor_mat)
plot(hclust(pheno_dist))
