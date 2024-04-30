#==============================
# get list of phenotypes with
# sample sizes >5000
# for genetic correlation
#==============================

library(data.table)
library(magrittr)

pheno <- fread("data/plink/adc_np.pheno", na.strings = "-1")

n_gt_5000 <- character(0)
for (i in colnames(pheno)[3:length(colnames(pheno))]) {
  if (pheno[!is.na(get(i)), .N] > 1400) {
    n_gt_5000 <- c(n_gt_5000, i)
  }
}

dt <- data.table(n_gt_5000)
fwrite(dt, file = "data/tmp/pheno_over_5000.txt", quote = FALSE, col.names = FALSE, sep = " ")

pheno_over_5000 <- pheno[, ..n_gt_5000]
pheno_over_5000 <- pheno_over_5000[complete.cases(pheno_over_5000)]

n_cols <- ncol(pheno_over_5000)
cor_mat <- matrix(nrow = n_cols, ncol = n_cols)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(pheno_over_5000)
for (i in 1L:n_cols) {
  pheno1 <- colnames(pheno_over_5000)[i]
  for (j in i:n_cols) {
    pheno2 <- colnames(pheno_over_5000)[j]
    cor_mat[i, j] <- cor_mat[j, i] <- pheno_over_5000[, cor(get(..pheno1), get(..pheno2))]
  }
}

cor_mat <- round(cor_mat, 2)

pheno1 <- pheno[, ..n_gt_5000]
cors <- cor(pheno1, use = "pairwise.complete.obs")
