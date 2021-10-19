#!/bin/Rscript --vanilla

#=====================================
# Format plink .assoc.logistic file
# into LDAK summary statistics format
# Predictor A1 A2 Z
#=====================================

pacman::p_load(data.table, magrittr, stringi, rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
bim_path <- args[1]
results_path <- args[2]
out_path <- args[3]
if (exists(args[4])) {
  if (args[4] %in% c("hg19", "hg38")) {
    build <- args[4]
  } else {
    stop("Build must be either: hg19, hg38")
  }
} else {
  build <- "hg38"
}

bim <- fread(bim_path)
results <- fread(results_path)

# # args while testing
# bim <- fread("data/plink/adc_np.bim")
# results <- fread("output/gwas_results/adc_np_HS.assoc.logistic")
# out_path <- "data/sumstats/HS_sumstats.txt"

# merge bim and results file to create sumstats data.table
results <- results[complete.cases(results)]
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
sumstats <- results[bim, on = c("CHR", "SNP", "BP", "A1")]
sumstats <- sumstats[complete.cases(sumstats)]

# LDAK only accepts SNPs for genetic correlation
snp_alleles <- c("A", "C", "T", "G")
sumstats <- sumstats[A1 %in% snp_alleles & A2 %in% snp_alleles]
sumstats[, ChrBP := paste0(CHR, ":", BP)]
sumstats <- sumstats[!(duplicated(ChrBP))]
sumstats[, seq_name := paste0("chr", CHR)]

# convert to hg19 build if current build is hg38
if (build == "hg38") {
  setnames(sumstats, "BP", "BP_hg38")
  
  # commented out because of issues installing GEOquery
  # copy over chain file
  if (!file.exists("raw_data/hg38ToHg19.over.chain")) {
    file.copy("/data_global/liftOver/hg38ToHg19.over.chain.gz",
              "raw_data/hg38ToHg19.over.chain.gz")
    GEOquery::gunzip("raw_data/hg38ToHg19.over.chain.gz",
                     "raw_data/hg38ToHg19.over.chain")
  }
  # convert to GRanges class
  seq.ranges <- IRanges(start = sumstats$BP_hg38,
                        end = sumstats$BP_hg38,
                        names = sumstats$SNP)
  sumstats.gr <- GRanges(
    seqnames = factor(sumstats$seq_name),
    ranges = seq.ranges,
    strand = factor(rep("*", nrow(sumstats)))
  )
  
  hg38.to.hg19.chain <- import.chain("raw_data/hg38ToHg19.over.chain")
  
  sumstats.gr.liftover <- liftOver(sumstats.gr, hg38.to.hg19.chain)
  sumstats.gr.liftover<- unlist(sumstats.gr.liftover)
  SNP <- names(sumstats.gr.liftover)
  sumstats.gr.liftover <- as.data.table(sumstats.gr.liftover)
  sumstats.gr.liftover[, SNP := ..SNP]
  setnames(sumstats.gr.liftover, "start", "BP_hg19")
  
  sumstats <- merge(sumstats, sumstats.gr.liftover[, .(SNP, BP_hg19)], "SNP")
  sumstats[, Predictor := paste0(CHR, ":", BP_hg19)]
} else {
  setnames(sumstats, "ChrBP", "Predictor")
}

setnames(sumstats, c("STAT", "NMISS"), c("Z", "n"))
sumstats <- sumstats[, .(Predictor, A1, A2, n, Z)]
sumstats <- sumstats[!(Predictor %in% Predictor[duplicated(Predictor)])]

fwrite(sumstats, file = out_path, quote = F, sep = " ")

rm(list = ls())
pacman::p_unload(all)
