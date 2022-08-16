library(data.table)
library(magrittr)
`%!in%` <- function(x, v) {!(x %in% v)}

pheno_ord <- fread("data/mega/mega_np_ord.pheno", na.strings = '-1')
pheno_bin <- fread("data/mega/mega_np.pheno")
phenotypes_ord <- colnames(pheno_ord)[3:ncol(pheno_ord)]
phenotypes_bin <- colnames(pheno_bin)[3:ncol(pheno_bin)]
phenotypes <- c(phenotypes_ord, phenotypes_bin)

results_files_ord <- paste0("output/gwas/mega/polmm/", 
                            colnames(pheno_ord)[3:ncol(pheno_ord)],
                            "_polmm_results.txt")
results_files_bin <- paste0("output/gwas/mega/saige/",
                            colnames(pheno_bin)[3:ncol(pheno_bin)],
                            "_saige_results.txt")
clump_files_ord <- paste0("output/gwas/mega/polmm/", 
                          colnames(pheno_ord)[3:ncol(pheno_ord)],
                          "_polmm_results.clumped")
clump_files_bin <- paste0("output/gwas/mega/saige/",
                          colnames(pheno_bin)[3:ncol(pheno_bin)],
                          "_saige_results.clumped")

get_results_list <- function(files) {
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  return(file_list)
}

# 1. collate binary and ordinal outcome results
results_ord <- get_results_list(results_files_ord)
for (i in 1:(ncol(pheno_ord) - 2)) {
  results_ord[[i]][, Phenotype := colnames(pheno_ord)[i + 2]]
}
names(results_ord) <- phenotypes_ord

results_bin <- get_results_list(results_files_bin)
for (i in 1:(ncol(pheno_bin) - 2)) {
  results_bin[[i]][, Phenotype := colnames(pheno_bin)[i + 2]]
  setnames(results_bin[[i]], 
           c("CHR", "MarkerID", "BETA", "p.value"),
           c("chr", "SNPID", "beta", "pval.spa"))
}

names(results_bin) <- phenotypes_bin

results <- c(results_ord, results_bin)
rm(results_ord, results_bin)

# 2. Create subsets of GWAS sumstats with P-values below threshold
threshold <- 1e-5
for (phenotype in phenotypes) {
assign(paste0(phenotype, "_top"),
       results[[phenotype]][chr != 19][pval.spa < threshold])
}

# 3. Merge subset sumstats for phenotype pairs where OR are in same direction
# and print the number of SNPs for each pair
for (i in 1:(length(phenotypes) - 1)) {
  pheno1 <- phenotypes[i]
  for (j in (i+1):length(phenotypes)) {
    pheno2 <- phenotypes[j]
    assign(paste0(pheno1, "_", pheno2),
           merge(get(paste0(pheno1, "_top")
           )[, c("chr", "SNPID", "beta", "pval.spa")], 
           get(paste0(pheno2, "_top")
               )[, c("chr", "SNPID", "beta", "pval.spa")], 
                 c('chr', 'SNPID'),
                 suffixes = c(paste0("_", pheno1), paste0("_", pheno2))))
    # assign(paste0(pheno1, "_", pheno2),
    #        get(paste0(pheno1, "_", pheno2))[
    #          (get(paste0('beta_', pheno1)) > 0 & 
    #             get(paste0('beta_', pheno2)) > 0) |
    #            (get(paste0('beta_', pheno1)) < 0 & 
    #               get(paste0('beta_', pheno2)) < 0)])
  }
}

for (i in 1:(length(phenotypes) - 1)) {
  pheno1 <- phenotypes[i]
  for (j in (i+1):length(phenotypes)) {
    pheno2 <- phenotypes[j]
    if (nrow(get(paste0(pheno1, "_", pheno2))[chr != 19]) > 0) {
      print(paste0(pheno1, "_", pheno2, " : ", 
                   nrow(get(paste0(pheno1, "_", pheno2))[chr != 19])))
    }}}
