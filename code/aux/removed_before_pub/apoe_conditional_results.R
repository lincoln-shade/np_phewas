library(data.table)
library(magrittr)
apoe_start <- 44905796
apoe_end <- 44909393
rad <- 5e5
caa_snp <- "rs4803778"
bim <- fread("data/mega/mega_np.bim") %>% 
  .[V1 == 19]
setnames(bim, c("V2", "V4"), c("SNPID", "pos"))
cerad_apoe <- fread("output/gwas/mega/polmm/cerad_apoe_polmm_results.txt",
                    check.names = T) %>% 
  merge(., bim[, .(SNPID, pos)], "SNPID")
cerad_apoe[pos > apoe_start - rad][pos < apoe_end + rad][
  pval.spa == min(pval.spa)
]
cerad_apoe[SNPID == "rs12721051"]
cerad_apoe[SNPID == caa_snp]

braak_apoe <- fread("output/gwas/mega/polmm/braak_apoe_polmm_results.txt",
                    check.names = T) %>% 
  merge(., bim[, .(SNPID, pos)], "SNPID")
braak_apoe[SNPID == "rs769449"]
braak_apoe[pos > apoe_start - rad][pos < apoe_end + rad][
  pval.spa == min(pval.spa)
]
braak_apoe[SNPID == caa_snp]

late_apoe <- fread("output/gwas/mega/polmm/late_apoe_polmm_results.txt",
                   check.names = T) %>% 
  merge(., bim[, .(SNPID, pos)], "SNPID")
late_apoe[pos > apoe_start - rad][pos < apoe_end + rad][
  pval.spa == min(pval.spa)
]
late_apoe[SNPID == "rs769449"]
late_apoe[SNPID == caa_snp]

caa_apoe <- fread("output/gwas/mega/polmm/caa_ord_apoe_polmm_results.txt",
                  check.names = T) %>% 
  merge(., bim[, .(SNPID, pos)], "SNPID")
caa_apoe[pos > apoe_start - rad][pos < apoe_end + rad][
  pval.spa == min(pval.spa)
]
caa_apoe[SNPID == "rs4420638"]
caa_apoe[SNPID == caa_snp]

abeta_apoe <- fread("output/gwas/mega/polmm/diffuse_abeta_apoe_polmm_results.txt",
                    check.names = T) %>% 
  merge(., bim[, .(SNPID, pos)], "SNPID")
abeta_apoe[pos > apoe_start - rad][pos < apoe_end + rad][
  pval.spa == min(pval.spa)
]
abeta_apoe[SNPID == caa_snp]