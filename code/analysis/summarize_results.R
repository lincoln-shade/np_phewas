

library(data.table)
library(magrittr)

file_name <- "output/act/ge_atherosclerosis_id_bin.assoc.logistic"
format_results <- function(file_name) {
  results <- fread(file_name)
  results <- results %>% 
    .[!(is.na(P))] %>% 
    setorder(P)
}

merge_results <- function(dt1, dt2, suf = c("_adc", "_act")) {
  merge(dt1, dt2, by = c("CHR", "BP", "A1", "SNP"), suffixes = suf)
}

