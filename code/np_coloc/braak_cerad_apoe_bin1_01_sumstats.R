
#=========================================================================
# perform colocalization analysis between Braak and CERAD at APOE and
# BIN1 loci
#=========================================================================

library(data.table)
library(magrittr)

pheno <- fread("data/mega/mega_np.pheno")
phenotypes <- c('braak56', 'cerad3')
braak_clumped <- fread("output/gwas/mega/braak56_clumped.clumped")
apoe_snp <- 'rs4420638'
bin1_snp <- 'rs6733839'
lead_snps <- c(apoe_snp, bin1_snp)
locus_radius_bp <- 2e5

file_path <- function(phenotype, prefix="output/gwas/mega") {
  paste0(prefix, '/', phenotype, '.assoc.logistic')
}

format_results <- function(file_name) {
  results <- fread(file_name)
  results <- results %>% 
    .[!(is.na(P))]
}

# 1. Read in GWAS sumstats for each phenotype
for (phenotype in phenotypes) {
  assign(phenotype, format_results(file_path(phenotype)))
}

# 2. subset sumstats to loci of interest
subset_locus <- function(dt, snp, radius=2e5) {
  snp_loc <- dt[SNP == snp][1, c(CHR, BP)]
  locus <- dt[CHR == snp_loc[1]
  ][
    (BP >= snp_loc[2] - radius) & (BP <= snp_loc[2] + radius)
    ]
}

braak56_apoe <- subset_locus(braak56, apoe_snp)
braak56_bin1 <- subset_locus(braak56, bin1_snp)
cerad3_apoe <- subset_locus(cerad3, apoe_snp)
cerad3_bin1 <- subset_locus(cerad3, bin1_snp)

coloc_sumstats <- list(braak56_apoe, cerad3_apoe, braak56_bin1, cerad3_bin1)

save(coloc_sumstats, file = "data/mega/coloc/apoe_bin1_coloc_sumstats.Rdata")
fwrite(braak56_apoe[, .(SNP)], 
       file = "data/mega/coloc/apoe_snps.txt",
       quote = FALSE, 
       col.names = FALSE,
       sep = ' ')
fwrite(braak56_bin1[, .(SNP)], 
       file = "data/mega/coloc/bin1_snps.txt",
       quote = FALSE, 
       col.names = FALSE,
       sep = ' ')
