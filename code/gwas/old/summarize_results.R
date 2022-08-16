
library(data.table)
library(magrittr)
source("code/functions/plink2_to_plink_format.R")

pheno <- fread("data/mega/mega_np.pheno", na.strings = '-1')
phenotypes <- colnames(pheno)[3:length(colnames(pheno))]

file_path <- function(phenotype, prefix="output/gwas/mega") {
  paste0(prefix, '/mega.', phenotype, '.glm.logistic')
}

format_results <- function(file_name) {
  results <- fread(file_name)
  results <- results %>% 
    .[!(is.na(P))] %>% 
    setorder(P)
}

for (phenotype in phenotypes) {
  assign(phenotype, format_results(file_path(phenotype)))
  assign(phenotype, plink2_to_plink_format(get(phenotype)))
  assign(paste0(phenotype, '_clumped'),
         fread(paste0("output/gwas/mega/mega.", phenotype, ".clumped")))
}

get_snp <- function(dt, chr, snp) {
  print(dt[CHR == chr][SNP == snp])
}

for (phenotype in phenotypes) {get_snp(get(phenotype), 2, 'rs6733839')}
for (phenotype in phenotypes) {print(phenotype); get_snp(get(phenotype), 7, 'rs4721058')}
for (phenotype in phenotypes) {print(phenotype); get_snp(get(phenotype), 19, 'rs4420638')}
