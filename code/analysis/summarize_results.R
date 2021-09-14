

library(data.table)
library(magrittr)

format_results <- function(phenotype) {
  results <- fread(paste0("output/gwas_results/adc_np_", phenotype, ".assoc.logistic"))
  results <- results %>% 
    .[!(is.na(P))] %>% 
    setorder(P)
}

pheno <- fread("data/plink/adc_np.pheno")

for (i in colnames(pheno)[3:ncol(pheno)]) {
  assign(i, format_results(i) %>% head())
}

for (i in colnames(pheno)[3:ncol(pheno)]) {
  if (nrow(get(i)) == 0) {
    rm(list = paste(i))
  }
}


