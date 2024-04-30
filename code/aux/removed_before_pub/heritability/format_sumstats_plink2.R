#!/usr/bin/Rscript --vanilla

#=====================================
# Format plink .assoc.logistic file
# into LDAK summary statistics format
# Predictor A1 A2 Z
#=====================================

library(data.table)
library(magrittr)
library(stringi)
library(rtracklayer)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--sumstats")
parser$add_argument("--N")
parser$add_argument("--out")
parser$add_argument("--build", default="hg38")
parser$add_argument("--tags", default="raw_data/ldak.thin.hapmap.gbr.tagging",
                    help="Tagging file for LDAK heritability model")
args <- parser$parse_args()

results_path <- args$sumstats
out_path <- args$out
tags_path <- args$tags

if (args$build %in% c("hg19", "hg38")) {
  build <- args$build
  print(paste("Variant positions are assumed to map to build:", build))
  if (build == "hg38") {
    print("Variants will be remapped to build hg19.")
  }
} else {
  stop("Build must be either: hg19, hg38")
}

results <- fread(results_path)
# results[Allele1 != ALT, T_STAT := -T_STAT]
results[, Z := Effect / StdErr]
results <- results[, .(`chr`, bp, MarkerName, Allele2, Allele1, OBS_CT, Z)]
setnames(results, colnames(results), c('CHR', 'BP', 'SNP', 'A2', 'A1', 'n', 'Z'))



# LDAK only accepts SNPs for genetic correlation
snp_alleles <- c("A", "C", "T", "G")
sumstats <- results[A1 %in% snp_alleles & A2 %in% snp_alleles]
sumstats[, ChrBP := paste0(CHR, ":", BP)]
sumstats <- sumstats[!(duplicated(ChrBP))]
sumstats[, seq_name := paste0("chr", CHR)]

# convert to hg19 build if current build is hg38
if (build == "hg38") {
  setnames(sumstats, "BP", "BP_hg38", skip_absent = TRUE)
  
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

sumstats <- sumstats[, .(Predictor, A1, A2, n, Z)]
sumstats <- sumstats[!(Predictor %in% Predictor[duplicated(Predictor)])]

tags <- fread(tags_path)
sumstats <- sumstats[Predictor %in% tags$Predictor]

fwrite(sumstats, file = out_path, quote = F, sep = " ")
