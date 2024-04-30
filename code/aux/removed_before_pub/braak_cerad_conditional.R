
library(data.table)
library(magrittr)
library(ordinal)
library(mediation)
library(MASS)

pheno <- setDT(readRDS('data/mega/mega_np.Rds'))[
  , .(FID, IID, braak, cerad, braak56, cerad3)
]

pheno[, braak_ord := as.ordered(braak)]
pheno[, cerad_ord := as.ordered(cerad)]
  
covar <- fread("data/mega/mega_np.covar")
snps <- fread("data/mega/conditional/bin1.raw")
related_rm <- fread("data/mega/related_rm/cerad3.remove", header = FALSE)
dt <- merge(pheno, covar, c('FID', 'IID'))
dt <- merge(dt, snps[, .(FID, IID, rs6733839_C)], c('FID', 'IID'))
dt <- dt[!(IID %in% related_rm$V2)]

dt[, age_death := scale(age_death)]
# dt[, age_death2 := scale(age_death2)]

# try a variety of models
covars <- paste(colnames(covar)[3:ncol(covar)], collapse = ' + ')

f_braak56 <- as.formula(paste('braak56 ~ rs6733839_C +', covars))
f_cerad3 <- as.formula(paste('cerad3 ~ rs6733839_C +', covars))
f1 <- as.formula(paste('braak56 ~ rs6733839_C + cerad3 +', covars))
f2 <- as.formula(paste('cerad3 ~ rs6733839_C + braak56 +', covars))
# f3 <- as.formula(paste('braak56 ~ rs4420638_G + cerad3 +', covars))
# f4 <- as.formula(paste('cerad3 ~ rs4420638_G + braak56 +', covars))

f5 <- as.formula(paste('braak56 ~ rs6733839_C + cerad +', covars))
f6 <- as.formula(paste('cerad3 ~ rs6733839_C + braak +', covars))
# f7 <- as.formula(paste('braak56 ~ rs4420638_G + cerad +', covars))
# f8 <- as.formula(paste('cerad3 ~ rs4420638_G + braak +', covars))


f9 <- as.formula(paste('braak_ord ~ rs6733839_C + cerad +', covars))
f10 <- as.formula(paste('cerad_ord ~ rs6733839_C + braak +', covars))
# f11 <- as.formula(paste('braak_ord ~ rs4420638_G + cerad +', covars))
# f12 <- as.formula(paste('cerad_ord ~ rs4420638_G + braak +', covars))
# 

f13 <- as.formula(paste('braak_ord ~ rs6733839_C + ', covars))
f14 <- as.formula(paste('cerad_ord ~ rs6733839_C + ', covars))


f15 <- as.formula(paste('braak_ord ~ rs6733839_C + cerad_ord +', covars))
f16 <- as.formula(paste('cerad_ord ~ rs6733839_C + braak_ord +', covars))

m1 <- dt[, glm(f1, family = binomial, data = .SD)]
summary(m1)

m2 <- dt[, glm(f2, family = binomial, data = .SD)]
summary(m2)

m3 <- dt[, glm(f3, family = binomial, data = .SD)]
summary(m3)

m4 <- dt[, glm(f4, family = binomial, data = .SD)]
summary(m4)

m_gwas <- dt[, glm(f_gwas, family = binomial, data = .SD)]
summary(m_gwas)

m5 <- dt[, glm(f5, family = binomial, data = .SD)]
summary(m5)

m6 <- dt[, glm(f6, family = binomial, data = .SD)]
summary(m6)

m7 <- dt[, glm(f7, family = binomial, data = .SD)]
summary(m7)

m8 <- dt[, glm(f8, family = binomial, data = .SD)]
summary(m8)

m9_null <- dt[, polr(f9_null, data = .SD, Hess = TRUE)]
summary(m9_null)
m9 <- dt[, polr(f9, data = .SD, Hess = TRUE)]
summary(m9)

m10_null <- dt[, polr(f10_null, data = .SD, Hess = TRUE)]
summary(m10_null)
m10 <- dt[, polr(f10, data = .SD, Hess = TRUE)]
summary(m10)


m13 <- dt[, clm(f13, data = .SD)]
summary(m13)
m14 <- dt[, clm(f14, data = .SD)]
summary(m14)

m15 <- dt[, clm(f15, data = .SD)]
summary(m15)


m16 <- dt[, clm(f16, data = .SD)]
summary(m16)

#----------------------------------------------
# Mediation analysis
#----------------------------------------------
# 1. Ordinal
keep_vars <- c('rs6733839_C', 'braak_ord', 'cerad_ord', colnames(covar)[3:ncol(covar)])
dt <- dt[, ..keep_vars][complete.cases(dt)]
iv_to_dv <- dt[, polr(f14, data = .SD, Hess = TRUE)]
iv_to_m <- dt[, polr(f13, data = .SD, Hess = TRUE)]
iv_m_to_dv <- dt[, polr(f16, data = .SD, Hess = TRUE)]

mr <- mediate(iv_to_m, 
              iv_m_to_dv, 
              mediator = 'braak_ord', 
              treat = 'rs6733839_C', 
              boot = TRUE)

save(mr, file = "output/gwas/mega/mediation_results.Rdata")


dt[braak_ord %in% as.character(0:4), braak56 := 0]
dt[braak_ord %in% as.character(5:6), braak56 := 1]
dt[cerad_ord %in% as.character(0:2), cerad3 := 0]
dt[cerad_ord %in% as.character(3), cerad3 := 1]

iv_to_dv <- dt[, glm(f_cerad3, data = dt, family = binomial())]
iv_to_m <- dt[, glm(f_braak56, data = dt, family = binomial())]
iv_m_to_dv <- dt[, glm(f2, data = dt, family = binomial())]

mr <- mediate(iv_to_m, 
              iv_m_to_dv, 
              mediator = 'braak56', 
              treat = 'rs6733839_C', 
              boot = TRUE)
