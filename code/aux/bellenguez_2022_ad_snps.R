library(data.table)
library(readxl)
library(stringi)

ad_snps_old <- read_xlsx('raw_data/bellenguez_et_al_2022_snps.xlsx')
ad_snps_new <- read_xlsx('raw_data/bellenguez_et_al_2022_snps.xlsx', sheet = 2)

braak <- fread('output/gwas/mega/braak_ord_results.txt')
braak[, SNP := stri_replace_first_regex(SNP, '_[ACTG]*', '')]
braak56 <- fread('output/gwas/mega/mega.braak56.glm.logistic')

npe_pc1 <- fread('tmp/mega.PC1.glm.linear')

cerad3 <- fread('output/gwas/mega/mega.cerad3.glm.logistic')
cerad <- fread('output/gwas/mega/braak_ord_results.txt')
braak[, SNP := stri_replace_first_regex(SNP, '_[ACTG]*', '')]

braak56_ad_snps_old <- braak56[ID %in% ad_snps_old$Variant]
braak[SNP %in% ad_snps_old$Variant]
braak[SNP %in% ad_snps_new$Variant]


npe_pc1_ad_snps_old <- npe_pc1[ID %in% ad_snps_old$Variant]
npe_pc1_ad_snps_new <- npe_pc1[ID %in% ad_snps_new$Variant]

cerad3_ad_snps_old <- cerad3[ID %in% ad_snps_old$Variant]
