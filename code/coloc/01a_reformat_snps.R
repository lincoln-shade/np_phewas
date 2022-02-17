#!/usr/bin/Rscript --vanilla

#=========================================
# take input variant ids and convert to 
# [chr]_[bp_hg38]_[ref]_[alt]_b38 format
#=========================================
options(warn=-1, message=-1)
library(data.table)
library(magrittr)
library(rentrez)
source("code/functions/strip_alleles.R")
variants <- fread("cat /dev/stdin", header = FALSE)$V1
args <- commandArgs(trailingOnly = TRUE)
phenotype <- args[1]
# variants <- fread(commandArgs(trailingOnly = T), header = F)$V1
if (is.null(variants)) {
  stop("No SNPs found!")
}
#-------------------------
# variants with rsIDs
#-------------------------
rsids <- variants[grep("rs", variants)]
# if (length(rsids) == 0) stop("You must input at least 1 rsID")

rsids <- stri_replace_first_regex(rsids, "rs", "") %>%
  strip_alleles() %>%
  as.integer()

# retrieve SNP summaries from dbSNP
dbsnp <- entrez_summary("snp", rsids, always_return_list = TRUE)

rsids <- as.character(rsids)

reformat_snp <- function(snp, dbsnp) {
  smry <- dbsnp[[snp]]
  alleles <- stri_extract_first_regex(smry[["docsum"]], "[:alpha:]*/[:alpha:]*") %>%
    gsub(pattern = "/", replacement = "_", x=.)
  pos <- gsub(":", "_", smry[["chrpos"]])
  snp.new.lab <- paste0("chr", pos, "_", alleles, "_b38")
  snp.new.lab
}

rsids_new <- rep(NA, length(rsids))
for (i in 1:length(rsids)) {
  rsids_new[i] <- reformat_snp(rsids[i], dbsnp)
}

#-----------------------------
# variants without rsIDs
#-----------------------------

no_rsid <- variants[grep("rs", variants, invert = TRUE)]
no_rsid_new <- paste0(stri_replace_all_fixed(no_rsid, ":", "_"), "_b38")

#-------------------------
# merge and write output
#------------------------

fwrite(as.data.table(c(rsids_new, no_rsid_new)), 
       file = paste0("data/tmp/snps_", phenotype, ".tmp"), 
       row.names = F, 
       col.names = F, 
       quote = F
)

rsid_key <- data.table("rsid"=c(rsids, no_rsid), 
                       "variant_id"=c(rsids_new, no_rsid_new)
)

fwrite(rsid_key, 
       file = paste0("data/tmp/rsid_key_", phenotype, ".tmp"), 
       row.names = F, 
       col.names = T, 
       quote = F
)

rsid_key[, cat_py := paste0(variant_id, "\n")]
if (nrow(rsid_key) > 0) {
  cat(rsid_key$cat_py, file = stdout())
} else {
  cat("No SNPs found!")
}

