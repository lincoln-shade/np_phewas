#=======================================================
# create data set with 4 CpG, SNP, CAA, and covariates
#=======================================================

library(data.table)
library(readxl)
library(ordinal)

load("data/cpg.RData")
covar <- setDT(read_xlsx("/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))
rosmap_clin <- fread("data/synapse/ROSMAP_clinical.csv")
rosmap_clin[, apoe := factor(apoe_genotype, 
                             levels = c(33, 34, 44, 23, 22), 
                             labels = c("apoe33", "apoe34", "apoe44", 
                                        "apoe23", "apoe22"))]
rosmap_clin[, projid := as.character(projid)]
# rs7247551 genotype
# cat 
# plink --bfile ../np_phewas/data/rosmap/rosmap_np --snp rs7247551 --recode A --out ../multi_omics/data/rs7247551
snp_geno <- fread("data/rs7247551.raw")
snp_geno[,  rs7247551_G := 2 - rs7247551_A]
snp_geno[, rs7247551_A := NULL]
rep_zeros <- function(x) {
    paste0(c(rep("0", 8 - nchar(x)), x), collapse = "")
}
rep_zeros <- Vectorize(rep_zeros)
snp_geno[, projid := rep_zeros(IID)]

cpg_apoe <- cpg[, .(IID, specimenID, conversion.efficiency, 
                    cg04401876, cg09555818, cg10169327, cg13119609)]
cpg_apoe[, IID := as.character(IID)]

rna <- fread("../np_phewas/data/rnaseq/rnaseq.txt.gz")
apoc2_col <- colnames(rna)[grep("ENSG00000234906", colnames(rna))]
apoe_col <- colnames(rna)[grep("ENSG00000130203", colnames(rna))]
rna_cols <- c("projid", "RIN", apoc2_col, apoe_col)
rna_apo <- rna[, ..rna_cols]
rna_apo[, projid := as.character(projid)]
setnames(rna_apo, c("ENSG00000234906.3", "ENSG00000130203.5"), c("APOC2", "APOE"))

# merge data into one data set
dt <- merge(covar, cpg_apoe, by.x = "projid", by.y = "IID", all.x = TRUE)
dt <- merge(dt, rosmap_clin[, .(projid, apoe)], "projid", all.x = TRUE)
dt <- merge(dt, snp_geno[, .(projid, rs7247551_G)], "projid", all.x = TRUE)
dt <- merge(dt, rna_apo, "projid", all.x = TRUE)
dt[, caa_4gp := as.ordered(caa_4gp)]
dt = dt[!is.na(caa_4gp)]

# write to .csv
fwrite(dt, file = "data/rosmap_data.csv", quote = FALSE)
