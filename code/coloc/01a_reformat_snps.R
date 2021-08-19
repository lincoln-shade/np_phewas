##---------------------------------------------
## take input rsIDs and convert to 
## [chr]_[bp_hg38]_[ref]_[alt]_b38 format
##---------------------------------------------

library(pacman)
p_load(data.table, magrittr, rentrez)
source("code/functions/strip_alleles.R")

rsids <- fread(commandArgs(trailingOnly = T), header = F)$V1
snps <- rsids
if (length(snps) == 0) stop("You must input at least 1 rsID")
snps <- stri_replace_first_regex(snps, "rs", "") %>% 
  strip_alleles() %>% 
  as.integer()

# retrieve SNP summaries from dbSNP
dbsnp <- entrez_summary("snp", snps)

snps <- as.character(snps)

reformat_snp <- function(snp, len) {
  if (length(snps) < 2) {smry <- dbsnp} else {smry <- dbsnp[[snps[i]]]}
  alleles <- stri_extract_first_regex(smry[["docsum"]], "[:alpha:]*/[:alpha:]*") %>% 
    gsub(pattern = "/", replacement = "_", x=.)
  pos <- gsub(":", "_", smry[["chrpos"]])
  snp.new.lab <- paste0("chr", pos, "_", alleles, "_b38")
  snp.new.lab
}

snps_new <- rep(NA, length(snps))
for (i in 1:length(snps)) {
  snps_new[i] <- reformat_snp(snps[i], length(snps))
}

fwrite(as.data.table(snps_new), file = "data/tmp/snps.tmp", row.names = F, col.names = F, quote = F)

rsid_key <- data.table("rsid"=rsids, "variant_id"=snps_new)

fwrite(rsid_key, file = "data/tmp/rsid_key.tmp", row.names = F, col.names = T, quote = F)

