library(data.table)
library(magrittr)
library(metaCCA)

hs <- fread("output/gwas_results/adc_np_HS.assoc.logistic")
late <- fread("output/gwas_results/adc_np_LATE.assoc.logistic")
frq <- fread("data/tmp/hs_late_snp.frq")
pheno <- fread("data/plink/adc_np.pheno", na.strings = "-1")
raw <- fread("data/tmp/hs_late_snp.raw")
# snp <- "rs7805419:T:C"
snp <- "rs254795:C:T"
hs <- hs[SNP == snp, .(SNP, A1, OR, SE)]
late <- late[SNP == snp, .(SNP, A1, OR, SE)]

hs[, hs_b := log(OR)]
late[, late_b := log(OR)]

sxy <- merge(hs, late, c("SNP", "A1"))
sxy[, allele_0 := "T"]
setnames(sxy, c("A1", "SE.x", "SE.y"), c("allele_1", "hs_se", "late_se"))
sxy <- data.frame(sxy[, .(allele_0, allele_1, hs_b, hs_se, late_b, late_se)], row.names = "rs1")
sxy$allele_0 <- factor(sxy$allele_0, levels = c("A", "C", "T", "G"))
sxy$allele_1 <- factor(sxy$allele_1, levels = c("A", "C", "T", "G"))

# elements of XY
hs_beta <- hs[SNP == snp, STAT / sqrt(NMISS)]
late_beta <- late[SNP == snp, STAT / sqrt(NMISS)]

# raw[, snp := scale(`rs7805419:T:C_C`)]

# covariance element in YY
pheno[!(is.na(HS)), hs_std := scale(HS)]
pheno[!(is.na(LATE)), late_std := scale(LATE)]
yy_cov <- pheno[, cov(hs_std, late_std, use = "p")]
# late_var <- pheno[, var(LATE, use = "c")]
# hs_var <- pheno[, var(HS, use = "c")]
# construct phenotype-phenotype covariance matrix
yy <- rbind(
  cbind(1, yy_cov),
  cbind(yy_cov, 1)
)

rownames(yy) <- colnames(yy) <- c("hs", "late")

# perform metaCCA
metaCcaGp(1, 
          S_XY = list(sxy),
          std_info = 0,
          S_YY = list(yy),
          N = 10000
)
