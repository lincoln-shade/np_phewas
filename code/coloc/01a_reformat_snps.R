#!/usr/bin/Rscript --vanilla

#=========================================
# take input variant ids and convert to 
# [chr]_[bp_hg38]_[ref]_[alt]_b38 format
#=========================================
options(warn=-1, message=-1)
library(data.table)
library(magrittr)
library(rentrez)
library(argparse)
source("code/functions/strip_alleles.R")

parser <- ArgumentParser()
parser$add_argument("--variant_file", "-f", 
                    help = "text file with SNPs (can also take from stdin)")
parser$add_argument("--out", "-o", help = "output files prefix")
args <- parser$parse_args()

if (!is.null(args$variant_file)) {
  variants <- fread(args$variant_file, header = FALSE) %>% 
    .[, V1]
} else {
  variants <- fread("cat /dev/stdin", header = FALSE) %>% 
    .[, V1]
}

if (length(variants) == 0) {
  stop("No variants input!")
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
max_var <- 200L # approx. # of SNPs one query can take w/o error
rsids_list <- split(rsids, ceiling(seq_along(rsids) / max_var))
dbsnp_list <- vector(mode = "list", length = length(rsids_list))
for (i in 1:length(dbsnp_list)) {
  dbsnp_list[[i]] <- entrez_summary("snp", 
                                  rsids_list[[i]], 
                                  always_return_list = TRUE)
}

dbsnp <- unlist(dbsnp_list, recursive = FALSE)
rsids <- as.character(rsids)

reformat_snp <- function(snp, dbsnp) {
  smry <- dbsnp[[snp]]
  alleles <- stri_extract_first_regex(smry[["docsum"]], 
                                      "[:alpha:]*/[:alpha:]*") %>%
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
if (length(no_rsid) > 0L) {
  no_rsid_new <- paste0(stri_replace_all_fixed(no_rsid, ":", "_"), "_b38")
} else {
  no_rsid_new <- character()
}


#-------------------------
# merge and write output
#------------------------

fwrite(as.data.table(c(rsids_new, no_rsid_new)), 
       file = paste0(args$out, "_snps.tmp"), 
       row.names = F, 
       col.names = F, 
       quote = F
)

rsid_key <- data.table("rsid"=c(rsids, no_rsid), 
                       "variant_id"=c(rsids_new, no_rsid_new)
)

fwrite(rsid_key, 
       file = paste0(args$out, "_key.tmp"), 
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

