#==================================================
# compare imputed APOE e diplotype vs synapse data
#==================================================

library(data.table)

rosmap_clin <- fread("data/synapse/ROSMAP_clinical.csv")
covar <- fread("~/projects/np_phewas/data/mega/mega_np_apoe.covar")

# harmonize variables
rosmap_clin[, apoe := paste0("apoe", apoe_genotype)]
rosmap_clin[, IID := as.character(projid)]

rosmap <- merge(covar, rosmap_clin[, .(IID, apoe_genotype)], "IID")

rosmap[, table(apoe, apoe_genotype)]

# does APOE genotype predict rs4803779 genotype
snp_geno <- fread("data/rs4803779.raw")
rep_zeros <- function(x) {
    paste0(c(rep("0", 8 - nchar(x)), x), collapse = "")
}
rep_zeros <- Vectorize(rep_zeros)
snp_geno[, IID := rep_zeros(IID)]

rosmap <- merge(rosmap, snp_geno[, .(IID, rs4803779_C)], "IID")
rosmap[, apoe_genotype := factor(apoe_genotype, levels = c(33, 34, 44, 23, 22))]
m_ld_null <- rosmap[!is.na(apoe_genotype), lm(rs4803779_C ~ 1)]
m_ld <- rosmap[, lm(rs4803779_C ~ apoe_genotype)]
summary(m_ld)
