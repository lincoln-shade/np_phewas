library(data.table)
library(ggplot2)
gwas = readRDS("output/coloc/meta/caa_apoe/GWAS/gwas_rs7247551.Rds")
cg13 = readRDS("output/coloc/meta/caa_apoe/ROSMAP/Brain_Prefrontal_Cortex_cg13119609_rs7247551.Rds")

cg13_pdata = merge(gwas$dataset, cg13$dataset, by.x = "rsID", by.y = "snp")
cg13_pdata[, ggplot(.SD, aes(-log10(pvalues.x), -log10(pvalues.y)))] +
  geom_point() +
  xlab("-log10(P_CAA)") +
  ylab("-log10(P_cg13")
