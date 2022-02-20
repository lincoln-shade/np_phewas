
library(data.table)
library(magrittr)
library(coloc)

pheno <- fread("data/mega/mega_np.pheno")
phenotypes <- c('braak56', 'cerad3')
load("data/mega/coloc/apoe_bin1_coloc_sumstats.Rdata")
braak56_apoe <- coloc_sumstats[[1]]
cerad3_apoe <- coloc_sumstats[[2]]
braak56_bin1 <- coloc_sumstats[[3]]
cerad3_bin1 <- coloc_sumstats[[4]]
rm(coloc_sumstats)
apoe_maf <- fread("data/mega/coloc/apoe.frq")
bin1_maf <- fread("data/mega/coloc/bin1.frq")

prop_cases <- function(phenotype) {
  n_cases <- pheno[get(phenotype) == 1, .N]
  n_controls <- pheno[get(phenotype) == 0, .N]
  proportion.cases <- n_cases / (n_cases + n_controls) 
}

braak56_prop_cases <- prop_cases('braak56')
cerad3_prop_cases <- prop_cases('cerad3')

p1 <- p2 <- 1e-4
p12 <-  1e-6
results_apoe <- coloc.abf(
  dataset1 = list(
    beta = log(braak56_apoe$OR),
    varbeta = (braak56_apoe$SE)^2,
    type = 'cc',
    s = braak56_prop_cases,
    N = braak56_apoe$NMISS,
    MAF = apoe_maf$MAF,
    snp = braak56_apoe$SNP
  ),
  dataset2 = list(
    beta = log(cerad3_apoe$OR),
    varbeta = (cerad3_apoe$SE)^2,
    type = 'cc',
    s = cerad3_prop_cases,
    N = cerad3_apoe$NMISS,
    MAF = apoe_maf$MAF,
    snp = cerad3_apoe$SNP
  ),
  p1 = p1,
  p2 = p2,
  p12 = p12
)

results_bin1 <- coloc.abf(
  dataset1 = list(
    beta = log(braak56_bin1$OR),
    varbeta = (braak56_bin1$SE)^2,
    type = 'cc',
    s = braak56_prop_cases,
    N = braak56_bin1$NMISS,
    MAF = bin1_maf$MAF,
    snp = braak56_bin1$SNP
  ),
  dataset2 = list(
    beta = log(cerad3_bin1$OR),
    varbeta = (cerad3_bin1$SE)^2,
    type = 'cc',
    s = cerad3_prop_cases,
    N = cerad3_bin1$NMISS,
    MAF = bin1_maf$MAF,
    snp = cerad3_bin1$SNP
  ),
  p1 = p1,
  p2 = p2,
  p12 = p12
)
