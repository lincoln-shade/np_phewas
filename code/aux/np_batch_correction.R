library(data.table)
library(magrittr)
library(psych)
library(ggplot2)

#-------------------------
# get complete cases
#-------------------------
bin_pheno <- fread("data/mega/mega_np.pheno", na.strings = "-1")
ord_pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = "-1")
pheno <- merge(ord_pheno, bin_pheno[, .(FID, IID, hs, grossinf, microinf)],
               by = c('FID', 'IID'))
covar <- fread("data/mega/mega_np.covar")

# remove TDP-43 since missing in ACT
# remove diffuse abeta since constructed differently in ROSMAP
pheno <- pheno[, -c("late", "diffuse_abeta")]

# keep complete cases since we want to compare distributions b/t data sources
pheno_complete <- pheno[complete.cases(pheno)]
pheno_complete_covar <- merge(pheno_complete, covar, c("FID", "IID"))

#-----------------------------------------
# perform linear models on all phenotypes
#-----------------------------------------
phenos <- colnames(pheno)[3:ncol(pheno)]
covars <- paste(colnames(covar)[3:ncol(covar)], collapse = ' + ')
pheno_residuals <- pheno_complete[, .(FID, IID)]

get_residuals <- function(npe) {
  f <- as.formula(paste(npe, "~", covars))
  m <- lm(f, data = pheno_complete_covar)
  res <- unname(m[["residuals"]])
  return(res)
}

for (i in phenos) {
  set(pheno_residuals,
      j = paste0(i, "_res"),
      value = get_residuals(i)
  )
}

#----------------------------------------------------------
# run PCA on residuals and compare dists. b/t data sources
#----------------------------------------------------------

pheno_pca <- pheno_residuals[, -c("FID", "IID")]

pca <- principal(pheno_pca, 
                 cor = "cor",
                 nfactors = ncol(pheno_pca),
                 rotate = "none") 

pca_scores <- pca$scores

pca_scores <- cbind(pheno_residuals[, .(FID, IID)], pca_scores)

# plot tops PCs and fill by FID
pc_plot <- function(pc, fid, dt, alpha=0.5) {
  pc_col <- sym(pc)
  fid_col <- sym(fid)
  ggplot(data = dt, aes(x = !!pc_col, fill=!!fid_col)) + 
    geom_density(alpha=alpha) +
    xlab(pc) + 
    theme_minimal() +
    scale_fill_discrete(labels=c("NACC", "ROSMAP", "ACT"))
}

for (i in 1:pca$factors) {
  assign(paste0("pc", i, "_plot"),
         pc_plot(paste0("PC", i), "FID", dt = pca_scores[FID != "ADNI"]))
  print(get(paste0("pc", i, "_plot")))
}

#--------------------------------------------------------------------
# run PCA on neuropath variables and compare dists. b/t data sources
#--------------------------------------------------------------------
pheno_pca_raw <- pheno_complete[, -c("FID", "IID")]

pca_raw <- principal(pheno_pca_raw, 
                 cor = "cor",
                 nfactors = ncol(pheno_pca_raw),
                 rotate = "none") 

pca_raw_scores <- pca_raw$scores

pca_raw_scores <- cbind(pheno_complete[, .(FID, IID)], pca_raw_scores)

for (i in 1:pca_raw$factors) {
  assign(paste0("pc", i, "_plot_raw"),
         pc_plot(paste0("PC", i), "FID", pca_raw_scores[FID != "ADNI"])
  )
  print(get(paste0("pc", i, "_plot_raw")))
}
