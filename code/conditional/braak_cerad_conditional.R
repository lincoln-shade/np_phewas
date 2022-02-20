
library(data.table)
library(readxl)
library(MASS)

pheno <- fread('data/mega/mega_np.pheno', na.strings = '-1')
covar <- fread("data/mega/mega_np.covar")
snps <- fread("data/mega/bin1_apoe_snps.raw")
related_rm <- fread("data/mega/related_rm/braak56.remove", header = FALSE)
dt <- merge(pheno, covar, c('FID', 'IID'))
dt <- merge(dt, snps[, .(FID, IID, rs6733839_T, rs4420638_G)], c('FID', 'IID'))
dt <- dt[!(IID %in% related_rm$V2)]

nacc <- setDT(readRDS('data/adc/np_qced.Rds'))
rosmap <- setDT(read_excel(
  "/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))
adni <- fread('raw_data/NEUROPATH_07_06_21.csv')
act <- setDT(read_xlsx("raw_data/act/DataForE235_20210421.xlsx",
                       sheet = 2))

# harmonize IDs between pheno file and individual study files
rosmap[, FID := '1']
setnames(rosmap, 'projid', 'IID')
act[, FID := '2']
setnames(act, 'IDfromDave', 'IID')
adni[, IID := as.character(RID)]
adni[, FID := 'ADNI']

varnames <- c('cerad', 'braak')
setnames(nacc, c('NACCNEUR', 'NACCBRAA'), varnames)
setnames(rosmap, c('braaksc'), c('braak'))
rosmap[ceradsc == 4, cerad := 0]
rosmap[ceradsc == 3, cerad := 1]
rosmap[ceradsc == 2, cerad := 2]
rosmap[ceradsc == 1, cerad := 3]
setnames(adni, c('NPNEUR', 'NPBRAAK'), varnames)
keep_vars <- c('FID', 'IID', 'braak', 'cerad')
covar_braak_cerad <- rbindlist(list(nacc[, ..keep_vars],
                                    rosmap[, ..keep_vars],
                                    act[, ..keep_vars],
                                    adni[, ..keep_vars]))

dt <- merge(dt, covar_braak_cerad, c('FID', 'IID'))
dt[, cerad_ord := ordered(cerad)]
dt[, braak_ord := ordered(braak)]

saveRDS(dt, file = "data/mega/conditional/dt.Rds")

# try a variety of models
covars <- paste(colnames(covar)[3:ncol(covar)], collapse = ' + ')
save(covars, file = "data/mega/conditional/covars.Rdata")
f_gwas <- as.formula(paste('braak56 ~ rs6733839_T +', covars))
f_gwas <- as.formula(paste('cerad3 ~ rs6733839_T +', covars))
f1 <- as.formula(paste('braak56 ~ rs6733839_T + cerad3 +', covars))
f2 <- as.formula(paste('cerad3 ~ rs6733839_T + braak56 +', covars))
f3 <- as.formula(paste('braak56 ~ rs4420638_G + cerad3 +', covars))
f4 <- as.formula(paste('cerad3 ~ rs4420638_G + braak56 +', covars))

f5 <- as.formula(paste('braak56 ~ rs6733839_T + cerad +', covars))
f6 <- as.formula(paste('cerad3 ~ rs6733839_T + braak +', covars))
f7 <- as.formula(paste('braak56 ~ rs4420638_G + cerad +', covars))
f8 <- as.formula(paste('cerad3 ~ rs4420638_G + braak +', covars))

f9_null <- as.formula(paste('braak_ord ~ cerad +', covars))
f9 <- as.formula(paste('braak_ord ~ rs6733839_T + cerad +', covars))
f10_null <- as.formula(paste('cerad_ord ~ braak +', covars))
f10 <- as.formula(paste('cerad_ord ~ rs6733839_T + braak +', covars))
f11 <- as.formula(paste('braak_ord ~ rs4420638_G + cerad +', covars))
f12 <- as.formula(paste('cerad_ord ~ rs4420638_G + braak +', covars))

f13_null <- as.formula(paste('braak_ord ~ ', covars))
f13 <- as.formula(paste('braak_ord ~ rs6733839_T + ', covars))
f14_null <- as.formula(paste('cerad_ord ~ ', covars))
f14 <- as.formula(paste('cerad_ord ~ rs6733839_T + ', covars))

f15_null <- as.formula(paste('braak_ord ~ cerad_ord +', covars))
f15 <- as.formula(paste('braak_ord ~ rs6733839_T + cerad_ord +', covars))
f16_null <- as.formula(paste('cerad_ord ~ braak_ord +', covars))
f16 <- as.formula(paste('cerad_ord ~ rs6733839_T + braak_ord +', covars))

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

m13_null <- dt[, polr(f13_null, data = .SD, Hess = TRUE)]
m13 <- dt[, polr(f13, data = .SD, Hess = TRUE)]
anova(m13_null, m13)

m14_null <- dt[, polr(f14_null, data = .SD, Hess = TRUE)]
m14 <- dt[, polr(f14, data = .SD, Hess = TRUE)]
anova(m14_null, m14)

m15_null <- dt[, polr(f15_null, data = .SD, Hess = TRUE)]
m15 <- dt[, polr(f15, data = .SD, Hess = TRUE)]
anova(m15_null, m15)
summary(m15)

m16_null <- dt[, polr(f16_null, data = .SD, Hess = TRUE)]
m16 <- dt[, polr(f16, data = .SD, Hess = TRUE)]
anova(m16_null, m16)
