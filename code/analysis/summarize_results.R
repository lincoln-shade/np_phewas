

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
    .[!(is.na(P))] %>% 
    setorder(P)
}

for (phenotype in phenotypes) {
  assign(phenotype, format_results(file_path(phenotype)))
}

for (phenotype in phenotypes) {print(head(get(phenotype)[CHR != 19]))}

get_snp <- function(dt, chr, snp) {
  dt[CHR == chr][SNP == snp]
}

for (phenotype in phenotypes) {get_snp(phenotype, 2, 'rs6733839')}
