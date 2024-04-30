library(data.table)
library(magrittr)
library(stringi)
rnaseq_t <- fread("data/mcbb/MayoRNAseq_RNAseq_TCX_geneCounts_normalized.tsv")
rnaseq <- t(rnaseq_t)
make_tidy <- function(m) {
  gene_id <- as.vector(m[1, ])
  specimenID <- rownames(m)[-1]
  m <- m[-1, ] %>% 
    apply(2, as.numeric) %>% 
    as.data.table()
  m[, specimenID := specimenID]
  setcolorder(m, "specimenID")
  
  colnames(m)[2:ncol(m)] <- gene_id
  m
}

rnaseq_tidy <- make_tidy(rnaseq)
rnaseq_tidy[, specimenID := stri_replace_last_fixed(specimenID, "_TCX", "")]

fwrite(rnaseq_tidy, file = "data/mcbb/rnaseq_gene_counts_normalized_tidy.csv",
       quote = FALSE)
