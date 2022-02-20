library(data.table)
library(magrittr)

pheno <- fread("data/mega/mega_np.pheno")
phenotypes <- colnames(pheno)[3:length(colnames(pheno))]

file_path <- function(phenotype, prefix="output/gwas/mega") {
  paste0(prefix, '/', phenotype, '.assoc.logistic')
}

format_results <- function(file_name) {
  results <- fread(file_name)
  results <- results %>% 
    .[!(is.na(P))]
}

# 1. Read in GWAS sumstats for each phenotype
for (phenotype in phenotypes) {
  assign(phenotype, format_results(file_path(phenotype)))
}

# 2. Create subsets of GWAS sumstats with P-values below threshold
threshold <- 1e-5
for (phenotype in phenotypes) {
assign(paste0(phenotype, "_top"),
       get(phenotype)[P < threshold])
}

# 3. Merge subset sumstats for phenotype pairs where OR are in same direction
# and print the number of SNPs for each pair
for (i in 1:(length(phenotypes) - 1)) {
  pheno1 <- phenotypes[i]
  for (j in (i+1):length(phenotypes)) {
    pheno2 <- phenotypes[j]
    assign(paste0(pheno1, "_", pheno2),
           merge(get(paste0(pheno1, "_top")), 
                 get(paste0(pheno2, "_top")), 
                 c('CHR', 'BP', 'SNP', 'A1', 'TEST'),
                 suffixes = c(paste0("_", pheno1), paste0("_", pheno2))))
    assign(paste0(pheno1, "_", pheno2),
           get(paste0(pheno1, "_", pheno2))[
             (get(paste0('OR_', pheno1)) > 1 & 
                get(paste0('OR_', pheno2)) > 1) |
               (get(paste0('OR_', pheno1)) < 1 & 
                  get(paste0('OR_', pheno2)) < 1)])
    print(paste0(pheno1, "_", pheno2, " : ", 
                 nrow(get(paste0(pheno1, "_", pheno2)))))
  }
}
