library(data.table)
library(ordinal)
library(MASS)
library(readxl)
rosmap_clin <- setDT(read_xlsx("/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))

###### CpG #####
load("data/cpg.RData")
cpg = cpg[, .(IID, specimenID, conversion.efficiency, cg09613507, cg23422036)]

np_cpg = merge(rosmap_clin, cpg, by.x = "projid", by.y = "IID", all.x = TRUE)
np_cpg[, late := as.ordered(tdp_st4)]

# LATE and HS mQTL models
m_late_cg09 = np_cpg[, clm(late ~ scale(cg09613507) + msex + scale(age_death) + scale(conversion.efficiency) + pmi, data = .SD)]
m_late_cg23 = np_cpg[, clm(late ~ scale(cg23422036) + msex + scale(age_death) + scale(conversion.efficiency) + pmi, data = .SD)]

m_hs_cg09 = np_cpg[, glm(hspath_typ ~ scale(cg09613507) + msex + scale(age_death) + scale(conversion.efficiency) + pmi, family = binomial)]
m_hs_cg23 = np_cpg[, glm(hspath_typ ~ scale(cg23422036) + msex + scale(age_death) + scale(conversion.efficiency) + pmi, family = binomial)]

summary(m_late_cg09)



##### RNA-Seq #####
rna <- fread("../np_phewas/data/rnaseq/rnaseq.txt.gz")
tmem106b_col <- colnames(rna)[grep("ENSG00000106460", colnames(rna))]
grn_col <- colnames(rna)[grep("ENSG00000030582", colnames(rna))]
rna_cols <- c("projid", "RIN", tmem106b_col, grn_col)
rna_ <- rna[, ..rna_cols]
rna_[, projid := as.character(projid)]
setnames(rna_, c(tmem106b_col, grn_col), c("TMEM106B", "GRN"))

np_cpg_rna = merge(np_cpg, rna_, "projid", all.x = TRUE)
late_snp = fread("data/rs2043539.raw")
late_snp[, IID := as.character(IID)]
np_cpg_rna = merge(np_cpg_rna, late_snp[, .(IID, rs2043539_A)], by.x = "projid", by.y = "IID", all.x = TRUE)
saveRDS(np_cpg_rna, file = "data/late_hs_qtl.Rds")

# HS <-> GRN
m_hs_grn = np_rna[, glm(hspath_typ ~ sqrt(GRN) + msex + age_death + pmi + RIN, family = binomial)]

# HS <-> TMEM106B
m_hs_tmem = np_rna[, glm(hspath_typ ~ sqrt(TMEM106B) + msex + age_death + pmi + RIN, family = binomial)]

# LATE <-> TMEM106B
m_late_tmem = np_cpg_rna[, clm(as.ordered(tdp_st4) ~ I(log(TMEM106B + 1)) + msex + age_death + pmi + RIN, data = .SD)]
m_tmem_late = np_rna[, lm(sqrt(TMEM106B) ~ as.factor(tdp_st4) + msex + age_death + pmi + RIN)]
##### CpG <-> TMEM106B #####
np_cpg_rna[, cor(sqrt(TMEM106B), cg09613507,use = "p")]

###### SNP analyses #####
###### plink --bfile ../np_phewas/data/rosmap/rosmap_np --snp rs2043539 --recode A --out ../multi_omics/data/rs2043539
###### plink --bfile ../np_phewas/data/rosmap/rosmap_np --snp rs7805419 --recode A --out ../multi_omics/data/rs7805419


hs_snp = fread("data/rs7805419.raw")
hs_snp[, IID := as.character(IID)]
np_cpg_rna = merge(np_cpg_rna, hs_snp[, .(IID, rs7805419_C)], by.x = "projid", by.y = "IID")

# late SNP 
m_snp_cg09 = np_cpg_rna[, lm(cg09613507 ~ rs2043539_A + msex + age_death + pmi + conversion.efficiency)]
m_snp_cg23 = np_cpg_rna[, lm(cg23422036 ~ rs2043539_A + msex + age_death + pmi + conversion.efficiency)]

m_snp_tmem = np_cpg_rna[, lm(sqrt(TMEM106B) ~ factor(rs2043539_A) + msex + age_death + pmi + RIN)]

# GRN snp
grn_snp = fread("data/rs5848.raw")
grn_snp[, IID := as.character(IID)]
np_rna_grn_snp = merge(np_rna, grn_snp[, .(IID, rs5848_T)], by.x = "projid", by.y = "IID")
m_snp_grn = np_rna_grn_snp[, lm(sqrt(GRN) ~ rs5848_T + msex + age_death)]
