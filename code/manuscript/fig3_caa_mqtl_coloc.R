library(data.table)
library(coloc)
library(ggplot2)
library(ggthemes)
caa_apoe <- fread("output/gwas/mega/polmm/caa_ord_apoe_polmm_results.txt")
pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = "-1")
covar <- fread("data/mega/mega_np_apoe.covar")
caa_apoe_snps <- caa_apoe[pval.spa < 5e-8, SNPID]
chr19_mqtl <- fread("/data_global/ROSMAP/xQTL_updated_data/sigmQTLs/sigmQTLchr19_50Kb.csv")
caa_apoe_mqtl <- chr19_mqtl[rsID %in% caa_apoe_snps]
chr19_all_mqtl <- fread("/data_global/ROSMAP/xQTL_updated_data/mQTLs/mQTLchr19_50Kb.csv")
caa_apoe_mqtl_loci <- chr19_all_mqtl[CpG %in% caa_apoe_mqtl$CpG]
caa_apoe <- caa_apoe[SNPID %in% caa_apoe_mqtl_loci$rsID]
caa_apoe_mqtl_loci <- caa_apoe_mqtl_loci[rsID %in% caa_apoe$SNPID]

setnames(caa_apoe, c("SNPID", "pval.spa"), c("rsID", "p"))
caa_apoe <- merge(caa_apoe, 
                  caa_apoe_mqtl_loci[!duplicated(rsID), 
                                     .(rsID, pos, A1, A2)], 
                  by = "rsID")
caa_apoe[, Phenotype := "CAA"]
setnames(caa_apoe_mqtl_loci, "CpG", "Phenotype")

# merge and plot
caa_mqtl_coloc <- rbind(caa_apoe[, .(Phenotype, rsID, pos, A1, A2, beta, p)],
                        caa_apoe_mqtl_loci[, .(Phenotype, rsID, pos, A1, A2, beta, p)])

ggp <- caa_mqtl_coloc[Phenotype %in% names(table(Phenotype))[c(1, 3:6)],
               ggplot(.SD, aes(pos / 1e6, -log10(p), color=factor(Phenotype)))] +
  geom_point() +
  theme_minimal() +
  xlab("Chromosome 19 Position (MB)") +
  ylab("-log10(p-value)") +
  labs(color = "Phenotype") +
  ggtitle("Colocalization of CAA and multiple mQTL in ROSMAP (PPH4 = 97%)") +
  theme(text = element_text(size = 16)) +
  scale_color_colorblind()

ggsave(plot = ggp,
       file = "doc/fig3_caa_coloc.png",
       units = "in",
       height = 8,
       width = 10)
  

# format for colocalization analysis
s <- pheno[caa_ord != 0, .N] / pheno[, .N]
caa_apoe[, s := ..s]
caa_apoe[, type := "cc"]
caa_apoe[, N := nrow(pheno)]
setnames(caa_apoe, c("p", "rsID"), c("pvalues", "snp"))
caa_apoe_mqtl_loci <- merge(caa_apoe_mqtl_loci, 
                            caa_apoe[, .(snp, MAF)],
                            by.x = "rsID",
                            by.y = "snp")
caa_apoe[, varbeta := (beta / qnorm(pval.norm/2, lower.tail = TRUE))^2]

dataset1 <- list(
  MAF=caa_apoe$MAF,
  type="cc",
  beta=caa_apoe$beta,
  varbeta=caa_apoe$varbeta,
  pvalues=caa_apoe$pvalues,
  snp=caa_apoe$snp,
  N = pheno[!is.na(caa_ord), .N],
  s = s
)
check_dataset(dataset1)

cpg1 <- caa_apoe_mqtl_loci[Phenotype == names(table(Phenotype))[5]]
dataset2 = list(
  beta = cpg1$beta,
  varbeta = cpg1$se^2,
  type = "quant",
  N = 468, # N for ROSMAP mQTL study
  MAF = cpg1$MAF,
  snp = cpg1$rsID
)
check_dataset(dataset2)

coloc.abf(dataset1, dataset2, p12=1e-5)
