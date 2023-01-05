library(data.table)
library(ggplot2)
library(forestplot)

make_95CI <- function(or, p, ci=0.95) {
  beta <- log(or)
  lt_zero <- beta < 0
  z <- qnorm(p/2, lower.tail = TRUE)
  sd <- beta / z
  zrad <- numeric(length = length(p))
  for (i in 1:length(p)){
  zrad[i] <- qnorm((1 - ci) / 2, lower.tail = lt_zero[i])
  }
  lower <- round(exp(beta + zrad * sd), 2)
  upper <- round(exp(beta - zrad * sd), 2)
  return(list(lower, upper))
}

#-----------------------------
# Braak stage PIK3R5 locus
#-----------------------------

braak_mega <- fread("output/gwas/mega/polmm/braak_polmm_results.txt")
braak_nacc <- fread("output/gwas/adc/polmm/braak_polmm_results.txt")
braak_rosmap <- fread("output/gwas/rosmap/polmm/braak_polmm_results.txt")
braak_act <- fread("output/gwas/act/polmm/braak_polmm_results.txt")
braak_adni <- fread("output/gwas/adni/polmm/braak_polmm_results.txt")
pik3r5_snp <- "rs72807981"
dt <- rbind(braak_nacc[chr == 17][SNPID == pik3r5_snp], 
                              braak_rosmap[chr == 17][SNPID == pik3r5_snp], 
                              braak_act[chr == 17][SNPID == pik3r5_snp], 
                              braak_adni[chr == 17][SNPID == pik3r5_snp],
                              braak_mega[chr == 17][SNPID == pik3r5_snp])
dt[, Study := c("NACC", "ROSMAP", "ACT", "ADNI", "Pooled")]
dt[, Lower := make_95CI(exp(beta), pval.spa)[[1]]]
dt[, Upper := make_95CI(exp(beta), pval.spa)[[2]]]
dt <- dt[Study != "ADNI"]
dt[, index := .N:1]

# 331x940 px
dt[, ggplot(data = .SD, aes(exp(beta), index, xmin=Lower, xmax=Upper))] +
  geom_point() +
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  geom_errorbarh(height=0.15) +
  geom_vline(xintercept = 1, linetype='dashed', alpha=.5) +
  scale_y_continuous(breaks = dt$index, labels = dt$Study) +
  ylab("Study") +
  xlab("Odds Ratio with 95% CI") # +
  # ggtitle("Braak stage results for PIK3R5 variant rs72807981, by study")

#-------------------------------
# Atherosclerosis COL4A1 locus
#-------------------------------
athero_mega <- fread("output/gwas/mega/polmm/atheroscler_polmm_results.txt")
athero_nacc <- fread("output/gwas/adc/polmm/atheroscler_polmm_results.txt")
athero_rosmap <- fread("output/gwas/rosmap/polmm/atheroscler_polmm_results.txt")
athero_act <- fread("output/gwas/act/polmm/atheroscler_polmm_results.txt")
athero_adni <- fread("output/gwas/adni/polmm/atheroscler_polmm_results.txt")
col4a1_snp <- "rs2000660"
dt <- rbind(athero_nacc[chr == 13][SNPID == col4a1_snp], 
            athero_rosmap[chr == 13][SNPID == col4a1_snp], 
            athero_act[chr == 13][SNPID == col4a1_snp], 
            athero_adni[chr == 13][SNPID == col4a1_snp],
            athero_mega[chr == 13][SNPID == col4a1_snp])
dt[, Study := c("NACC", "ROSMAP", "ACT", "ADNI", "Pooled")]
dt[, Lower := make_95CI(exp(beta), pval.spa)[[1]]]
dt[, Upper := make_95CI(exp(beta), pval.spa)[[2]]]
dt <- dt[Study != "ADNI"]
dt[, index := .N:1]

dt[, ggplot(data = .SD, aes(exp(beta), index, xmin=Lower, xmax=Upper))] +
  geom_point() +
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  geom_errorbarh(height=0.1) +
  geom_vline(xintercept = 1, linetype='dashed', alpha=.5) +
  scale_y_continuous(breaks = dt$index, labels = dt$Study) +
  ylab("Study") +
  xlab("Odds Ratio with 95% CI") # +
  # xlim(0, 5) +
  # ggtitle("Atherosclerosis results for COL4A1 variant rs2000660, by study")

#-------------------------------
# CAA APOC2 locus
#-------------------------------

caa_mega <- fread("output/gwas/mega/polmm/caa_ord_polmm_results.txt")
caa_nacc <- fread("output/gwas/adc/polmm/caa_ord_polmm_results.txt")
caa_rosmap <- fread("output/gwas/rosmap/polmm/caa_ord_polmm_results.txt")
caa_act <- fread("output/gwas/act/polmm/caa_ord_polmm_results.txt")
caa_adni <- fread("output/gwas/adni/polmm/caa_ord_polmm_results.txt")
caa_snp <- "rs4803774"
dt <- rbind(caa_nacc[chr == 19][SNPID == caa_snp], 
            caa_rosmap[chr == 19][SNPID == caa_snp], 
            caa_act[chr == 19][SNPID == caa_snp], 
            caa_adni[chr == 19][SNPID == caa_snp],
            caa_mega[chr == 19][SNPID == caa_snp])
dt[, Study := c("NACC", "ROSMAP", "ACT", "ADNI", "Pooled")]
dt[, Lower := make_95CI(exp(beta), pval.spa)[[1]]]
dt[, Upper := make_95CI(exp(beta), pval.spa)[[2]]]
dt <- dt[Study != "ADNI"]
dt[, index := .N:1]

dt[, ggplot(data = .SD, aes(exp(beta), index, xmin=Lower, xmax=Upper))] +
  geom_point() +
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  geom_errorbarh(height=0.1) +
  geom_vline(xintercept = 1, linetype='dashed', alpha=.5) +
  scale_y_continuous(breaks = dt$index, labels = dt$Study) +
  ylab("Study") +
  xlab("Odds Ratio with 95% CI") +
# xlim(0, 5) +
  ggtitle("CAA results for APOC2 variant rs4803774, by study")
