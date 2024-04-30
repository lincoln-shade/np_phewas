
#=========================================================================
# perform colocalization analysis between Braak and CERAD BIN1 locus
#=========================================================================

library(data.table)
library(magrittr)
library(coloc)

pheno <- fread("shared/nacc_rosmap_act_np_dichotomized.pheno")
phenotypes <- c('braak', 'cerad')
braak_clumped <- fread("output/gwas/metal/results/braak1_clump.clumped")
bin1_snp <- 'rs6733839'
locus_radius_bp <- 2e5

# 1. Read in GWAS sumstats for each phenotype
braak = fread("output/gwas/metal/results/braak1.csv")
cerad = fread("output/gwas/metal/results/cerad1.csv")

# 2. subset sumstats to loci of interest
subset_locus <- function(dt, snp, radius=2e5) {
  snp_loc <- dt[MarkerName == snp][1, c(chr, bp)]
  locus <- dt[chr == snp_loc[1]
  ][
    (bp >= snp_loc[2] - radius) & (bp <= snp_loc[2] + radius)
    ]
}

braak_bin1 <- subset_locus(braak, bin1_snp)
cerad_bin1 <- subset_locus(cerad, bin1_snp)

# 3. Get proportion of cases for each phenotype

prop_cases <- function(phenotype) {
  n_cases <- pheno[get(phenotype) == 1, .N]
  n_controls <- pheno[get(phenotype) == 0, .N]
  proportion.cases <- n_cases / (n_cases + n_controls) 
}

braak_prop_cases <- prop_cases('braak')
cerad_prop_cases <- prop_cases('cerad')

# 4. set coloc priors

p1 <- p2 <- 1e-4
p12 <-  1e-5

# 5. perform colocalization

results_bin1 <- coloc.abf(
  dataset1 = list(
    beta = braak_bin1$Effect,
    varbeta = (braak_bin1$StdErr)^2,
    type = 'cc',
    s = braak_prop_cases,
    N = pheno[!is.na(braak), .N],
    MAF = braak_bin1$Freq1,
    snp = braak_bin1$MarkerName
  ),
  dataset2 = list(
    beta = cerad_bin1$Effect,
    varbeta = (cerad_bin1$StdErr)^2,
    type = 'cc',
    s = cerad_prop_cases,
    N = pheno[!is.na(cerad), .N],
    MAF = cerad_bin1$Freq1,
    snp = cerad_bin1$MarkerName
  ),
  p1 = p1,
  p2 = p2,
  p12 = p12
)

# Results
# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 2.18e-07  1.95e-04  1.84e-06  6.51e-04  9.99e-01 
# [1] "PP abf for shared variant: 99.9%"