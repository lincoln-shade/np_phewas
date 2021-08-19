##=======================================
## prepare b-asc data for coloc analysis
##=======================================

library(pacman)
p_load(data.table, magrittr, rtracklayer, GEOquery, GENESIS, GWASTools, qusage, SeqArray, SNPRelate)
source("code/functions/strip_alleles.R")

# commandline arguments
# cargs 1 = phenotype file (remember phenotype needs to be binary  (1, 2) and third column)
# cargs 2 = minor allele frequency file .frq
# cargs 3 = plink regression file .assoc.logistic
# cargs 4 = chromosome integer
cargs <- commandArgs(trailingOnly = T)

##------------------
## B-ASC GWAS data
##------------------

# GWAS pheno file
pheno <- fread(cargs[1])
# proportion of participants who were cases
proportion.cases <- pheno[pheno[[3]] == 2, .N] / pheno[, .N]

# MAF data
gwas.maf <- fread(paste0(cargs[2]))
gwas.maf <- gwas.maf[, .(SNP, MAF, A1)]
gwas.maf[, SNP := strip_alleles(SNP)]

# GWAS summary stats
gwas <- fread(cargs[3])
gwas <- gwas[CHR == as.integer(cargs[4])]
setnames(gwas, "BP", "BP_hg19")
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

##---------------------------------------------------
## Convert variant hg19 positions to hg38 positions
##---------------------------------------------------

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

# keep only required columns
gwas <- gwas[, .(SNP, A1, CHR, BP_hg38, BP_hg19, MAF, P, beta, var.beta, NMISS, type, s)]
setnames(gwas, c("MAF", "NMISS"), c("maf", "n"))
gwas <- unique(gwas)


save(gwas, file = paste0("data/tmp/chr", cargs[4], "_gwas_sumstats.RData"))

# rm(list=ls())
# p_unload(all)
