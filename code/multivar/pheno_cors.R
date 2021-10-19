#==============================
# get list of phenotypes with
# sample sizes >5000
# for genetic correlation
#==============================

library(data.table)
library(magrittr)

bin_pheno <- fread("data/plink/adc_np.pheno", na.strings = "-1")
ord_pheno <- fread("data/adc_np_ord.txt", na.strings = "-1")

for (i in colnames(bin_pheno[, 3:ncol(bin_pheno)])) {
  if (is.numeric(bin_pheno[[i]])) {
    if (bin_pheno[!(is.na(get(i))), 
                  sum(get(..i))] %in% c(0, .N)) {
      bin_pheno[, (i) := NULL]
    }
  }
}

set_null <- c("NACCPICK", "NACCCBD", "NACCPROG", "NPLINF", "NPLAC", "NPHEM", "NPMICRO", "NPART", "NPOANG", "NPSCL")
for (i in set_null) {
  bin_pheno[, (i) := NULL]
}

bin_pheno_4k <- character()
for (i in colnames(bin_pheno)) {
  if (bin_pheno[!(is.na(get(i))), .N] >= 4000) {
    bin_pheno_4k <- c(bin_pheno_4k, i)
  }
}
bin_pheno_4k <- bin_pheno[, ..bin_pheno_4k]
bin_pheno_4k <- bin_pheno_4k[complete.cases(bin_pheno_4k)]

n_cols <- ncol(bin_pheno_4k) - 2L
cor_mat <- matrix(nrow = n_cols, ncol = n_cols)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(bin_pheno_4k[, -c("FID", "IID")])
for (i in 1L:n_cols) {
  pheno1 <- colnames(bin_pheno_4k)[i + 2L]
  for (j in i:n_cols) {
    pheno2 <- colnames(bin_pheno_4k)[j + 2L]
    cor_mat[i, j] <- cor_mat[j, i] <- bin_pheno_4k[, cor(get(..pheno1), get(..pheno2), use = "c")]
  }
}

cor_mat[is.na(cor_mat)] <- 0
cor_mat <- round(cor_mat, 2)

pca1 <- prcomp(bin_pheno_4k[, -c("FID", "IID")], scale. = T, center = T)
pcs <- data.table(FID = bin_pheno_4k$FID, 
                  IID = bin_pheno_4k$IID, 
                  pca1$x)

bin_pheno_4k[IID %in% pcs[PC1 > 4, IID]]
