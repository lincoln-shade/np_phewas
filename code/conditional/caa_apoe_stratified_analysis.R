# APOE-stratified analysis for CAA

library(data.table)
library(ordinal)

###### NACC ######
pheno <- fread("data/adc/adc_np_ord.txt", na.strings = "-1")
covar <- fread("data/adc/adc_np_apoe.covar")
# plink --bfile data/adc/adc_np --recode A --snp rs7247551 --out data/adc/apoe_snp_rs7247551
snp <- fread("data/adc/apoe_snp_rs7247551.raw")

dt <- Reduce(
  merge, 
  list(pheno, 
       covar, 
       snp[, .(FID, IID, 'snp' = rs7247551_G)]
  )
)

dt[, apoe := factor(
  apoe, 
  levels = c("apoe33", "apoe34", "apoe44", "apoe23", "apoe22", "apoe24")
)]

# analysis
make_formula <- function(v, covars) {
  formula(paste(c(v, covars), collapse = " + "))
}

covars <- colnames(covar)[3:(ncol(covar) - 1)]
f1 <- make_formula("caa ~ snp", covars)
dt[, caa := as.ordered(caa)]
# apoe33 OR = 1.18, P = 0.00039
m_apoe33 <- dt[apoe == "apoe33", clm(f1, data = .SD)]

# apoe34 OR = 1.28, P = 1.63e-6
m_apoe34 <- dt[apoe == "apoe34", clm(f1, data = .SD)]

# apoe44 OR = 1.20, P = 0.077
m_apoe44 <- dt[apoe == "apoe44", clm(f1, data = .SD)]

# apoe23 OR = 1.46, P = 0.0073
m_apoe23 <- dt[apoe == "apoe23", clm(f1, data = .SD)]

# apoe24 OR = 1.50, P = 0.079
m_apoe24 <- dt[apoe == "apoe24", clm(f1, data = .SD)]

# combined model with interaction effect
# no snp x APOExx beta significant
f2 <- make_formula("caa ~ snp*apoe", covars)
m_apoe_xsnp <- dt[, clm(f2, data = .SD)]

f2_null <- make_formula("caa ~ snp + apoe", covars)
m_apoe_xsnp_null <- dt[!is.na(caa), clm(f2_null, data = .SD)]

# chisq P = 4.1e-9
anova(m_apoe_xsnp_null, m_apoe_xsnp)

###### ROSMAP ######
pheno_rosmap <- fread("data/rosmap/rosmap_np.pheno", na.strings = "-1")
covar_rosmap <- fread("data/rosmap/rosmap_np_apoe.covar")
# plink --bfile data/rosmap/rosmap_np --recode A --snp rs7247551 --out data/rosmap/apoe_snp_rs7247551
snp_rosmap <- fread("data/rosmap/apoe_snp_rs7247551.raw")
snp_rosmap[, rs7247551_G := 2 - rs7247551_A]
snp_rosmap[, rs7247551_A := NULL]

dt_rosmap <- Reduce(
  merge, 
  list(pheno_rosmap, 
       covar_rosmap, 
       snp_rosmap[, .(FID, IID, 'snp' = rs7247551_G)]
  )
)

dt_rosmap[, apoe := factor(
  apoe, 
  levels = c("apoe33", "apoe34", "apoe44", "apoe23", "apoe22", "apoe24")
)]

covars_rosmap <- colnames(covar_rosmap)[3:(ncol(covar_rosmap) - 1)]
f1_rosmap <- make_formula("caa ~ snp", covars_rosmap)
f1a_rosmap = make_formula("caa ~ snp + apoe", covars_rosmap)
dt_rosmap[, caa := as.ordered(caa)]
# apoe33 OR = 0.79, P = 0.021
m_apoe33_rosmap <- dt_rosmap[apoe == "apoe33", clm(f1_rosmap, data = .SD)]

# apoe34 or apoe44 or apoe24 (combined due to low sample size of apoe44 (n=20) and apoe24 (n=19)) OR = 0.68, P = 0.015
m_apoe34_44_24_rosmap <- dt_rosmap[apoe %in% c("apoe34", "apoe44", "apoe24"), clm(f1a_rosmap, data = .SD)]

# apoe23 OR = 0.65, P = 0.13
m_apoe23_rosmap <- dt_rosmap[apoe == "apoe23", clm(f1_rosmap, data = .SD)]

# combined model with interaction effect
# no snp x APOExx beta significant
f2_rosmap <- make_formula("caa ~ snp*apoe", covars_rosmap)
m_apoe_xsnp_rosmap <- dt_rosmap[, clm(f2_rosmap, data = .SD)]

f2_null_rosmap <- make_formula("caa ~ snp + apoe", covars_rosmap)
m_apoe_xsnp_null_rosmap <- dt_rosmap[!is.na(caa), clm(f2_null_rosmap, data = .SD)]

# chisq P = 0.13
anova(m_apoe_xsnp_null_rosmap, m_apoe_xsnp_rosmap)


###### ACT ######
pheno_act <- fread("data/act/act_np.pheno", na.strings = "-1")
covar_act <- fread("data/act/act_np_apoe.covar")
# plink --bfile data/act/act_np --recode A --snp rs7247551:G:A --out data/act/apoe_snp_rs7247551
snp_act <- fread("data/act/apoe_snp_rs7247551.raw", check.names = T)
snp_act[, rs7247551_G := 2 - rs7247551.G.A_A]
snp_act[, rs7247551.G.A_A := NULL]

dt_act <- Reduce(
  merge, 
  list(pheno_act, 
       covar_act, 
       snp_act[, .(FID, IID, 'snp' = rs7247551_G)]
  )
)

dt_act[, apoe := factor(
  apoe, 
  levels = c("apoe33", "apoe34", "apoe44", "apoe23", "apoe22", "apoe24")
)]

covars_act <- colnames(covar_act)[3:(ncol(covar_act) - 1)]
f1_act <- make_formula("caa ~ snp", covars_act)
f1a_act = make_formula("caa ~ snp + apoe", covars_act)
dt_act[, caa := as.ordered(caa)]
# apoe33 OR = 0.88, P = 0.41
m_apoe33_act <- dt_act[apoe == "apoe33", clm(f1_act, data = .SD)]

# apoe34 or apoe44 or apoe24 (combined due to low sample size of apoe44 (n=20) and apoe24 (n=19)) OR = 0.68, P = 0.015
m_apoe34_44_24_act <- dt_act[apoe %in% c("apoe34", "apoe44", "apoe24"), clm(f1a_act, data = .SD)]

# apoe23 OR = 0.93, P = 0.72
m_apoe23_act <- dt_act[apoe == "apoe23", clm(f1_act, data = .SD)]

# combined model with interaction effect
# no snp x APOExx beta significant
f2_act <- make_formula("caa ~ snp*apoe", covars_act)
m_apoe_xsnp_act <- dt_act[, clm(f2_act, data = .SD)]

f2_null_act <- make_formula("caa ~ snp + apoe", covars_act)
m_apoe_xsnp_null_act <- dt_act[!is.na(caa), clm(f2_null_act, data = .SD)]

# chisq P = 0.91
anova(m_apoe_xsnp_null_act, m_apoe_xsnp_act)
