#!/usr/bin/env Rscript --vanilla

##=======================================
## prepare b-asc data for coloc analysis
##=======================================

library(data.table)
library(magrittr)
library(rtracklayer)
library(GEOquery)
library(GENESIS)
library(GWASTools)
library(qusage)
library(SeqArray)
library(SNPRelate)
source("code/functions/strip_alleles.R")

# commandline arguments
# args 1 = phenotype file (phenotype needs to be binary (0/1) 
# args 2 = minor allele frequency file .frq
# args 3 = plink regression file .assoc.logistic
# args 4 = chromosome integer
# args 5 = gwas phenotype name
args <- commandArgs(trailingOnly = T)

##------------------
## B-ASC GWAS data
##------------------

# GWAS pheno file
pheno <- fread(args[1], na.strings = '-1')
# proportion of participants who were cases
gwas_phenotype <- args[5]
n_cases <- pheno[get(gwas_phenotype) == 1, .N]
n_controls <- pheno[get(gwas_phenotype) == 0, .N]
proportion.cases <- n_cases / (n_cases + n_controls) 

# MAF data
gwas.maf <- fread(paste0(args[2]))
gwas.maf <- gwas.maf[, .(SNP, MAF, A1)]
gwas.maf[, SNP := strip_alleles(SNP)]

# GWAS summary stats
gwas <- fread(args[3])
gwas <- gwas[CHR == as.integer(args[4])]
gwas[, seq.name := paste0("chr", CHR)]
gwas[, SNP := strip_alleles(SNP)]
gwas[, beta := log(OR)]
gwas[, var.beta := SE**2]
# B-ASC GWAS was case-control study
gwas[, type := "cc"]
# proportion of participants cases
gwas[, s := proportion.cases]

gwas <- merge(gwas, gwas.maf, c("SNP", "A1"))
gwas <- gwas[grep("rs", SNP)] # remove non-rsID variants
rm(gwas.maf, pheno)

#---------------------------------------------------
# Convert variant hg19 positions to hg38 positions
# if needed
#---------------------------------------------------

build <- 38
if (build == 19) {
  setnames(gwas, "BP", "BP_hg19")
  
  # commented out because of issues installing GEOquery
  # copy over chain file
  if (!file.exists("raw_data/hg19ToHg38.over.chain")) {
    file.copy("/data_global/liftOver/hg19ToHg38.over.chain.gz",
              "raw_data/hg19ToHg38.over.chain.gz")
    GEOquery::gunzip("/raw_data/hg19ToHg38.over.chain.gz",
                     "raw_data/hg19ToHg38.over.chain")
  }
  # convert to GRanges class
  seq.ranges <- IRanges(start = gwas$BP_hg19,
                        end = gwas$BP_hg19,
                        names = gwas$SNP)
  gwas.gr <- GRanges(
    seqnames = factor(gwas$seq.name),
    ranges = seq.ranges,
    strand = factor(rep("*", nrow(gwas)))
  )
  
  hg19.to.hg38.chain <- import.chain("raw_data/hg19ToHg38.over.chain")
  
  gwas.gr.liftover <- liftOver(gwas.gr, hg19.to.hg38.chain)
  gwas.gr.liftover<- unlist(gwas.gr.liftover)
  SNP <- names(gwas.gr.liftover)
  gwas.gr.liftover <- as.data.table(gwas.gr.liftover)
  gwas.gr.liftover[, SNP := ..SNP]
  setnames(gwas.gr.liftover, "start", "BP_hg38")
  
  gwas <- merge(gwas, gwas.gr.liftover[, .(SNP, BP_hg38)], "SNP")
  gwas <- gwas[, .(SNP, A1, CHR, BP_hg38, BP_hg19, MAF, P, beta, var.beta, 
                   NMISS, type, s)]
} else {
  setnames(gwas, "BP", "BP_hg38")
  gwas <- gwas[, .(SNP, A1, CHR, BP_hg38, MAF, P, beta, var.beta, NMISS, 
                   type, s)]
}

# keep only required columns
setnames(gwas, c("MAF", "NMISS"), c("maf", "n"))
gwas <- unique(gwas)


save(gwas, file = paste0("data/tmp/chr", args[4], "_gwas_sumstats_", 
                         gwas_phenotype, ".RData"))

