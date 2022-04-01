library(data.table)
library(magrittr)
library(MASS)
library(mediation)
library(ordinal)

pheno <- setDT(readRDS('data/mega/mega_np.Rds'))[
  , .(FID, IID, braak, braak56, cerad, cerad3)
]
covar <- fread('data/mega/mega_np.covar')
covar[, msex_age_death := msex * age_death]
related_rm <- fread('data/mega/related_rm/cerad3.remove')
bin1_snps <- fread('data/mega/conditional/bin1.raw')
bin1_top_braak56_snp <- 'rs6733839'
bin1_top_snp_col <- colnames(bin1_snps)[
  grep(bin1_top_braak56_snp, colnames(bin1_snps))]

covars <- colnames(covar)[3:ncol(covar)] # %>%
  # .[!(. %in% c('act', 'ACT2', 'ACT3', 'adni'))]

dt <- Reduce(merge, 
             list(pheno, 
                  covar, 
                  bin1_snps[, .(FID, IID, 'snp' = get(..bin1_top_snp_col))]))

dt <- dt[!(IID %in% related_rm$V2)] %>% 
  .[complete.cases(.)]

# binary mediation analysis
f1 <- formula(paste0(c('cerad3 ~ snp', covars), collapse = ' + '))
f2 <- formula(paste0(c('braak56 ~ snp', covars), collapse = ' + '))
f3 <- formula(paste0(c('cerad3 ~ snp', 'braak56', covars), collapse = ' + '))
iv_to_dv <- dt[, glm(f1, data = .SD, family = binomial())]
iv_to_m <- dt[, glm(f2, data = .SD, family = binomial())]
iv_m_to_dv <- dt[, glm(f3, data = .SD, family = binomial())]

mr_binary <- mediate(iv_to_m, 
                     iv_m_to_dv, 
                     mediator = 'braak56', 
                     treat = 'snp', 
                     boot = TRUE, 
                     sims = 100)

# ordinal late mediation analysis
dt[, braak_ord := ordered(braak)]
dt[, cerad_ord := ordered(cerad)]
f4 <- formula(paste0(c('braak_ord ~ snp', covars), collapse = ' + '))
f5 <- formula(paste0(c('cerad_ord ~ snp', 'braak_ord', covars), collapse = ' + '))
f6 <- formula(paste0(c('cerad_ord ~ snp', covars), collapse = ' + '))
iv_to_m <- dt[, polr(f4, data = .SD, Hess = TRUE)]
iv_m_to_dv <- dt[, polr(f5, data = .SD, Hess = TRUE)]
iv_to_dv <- dt[, polr(f6, data = .SD, Hess = TRUE)]

mr_ord <- mediate(iv_to_m, 
                  iv_m_to_dv, 
                  mediator = 'braak_ord', 
                  treat = 'snp', 
                  boot = TRUE, 
                  sims = 100)

# some more conditional analyses
m1 <- dt[braak < 4, clm(f6, data = .SD)]
summary(m1)

m2 <- dt[braak >= 4, clm(f6, data = .SD)]
summary(m2)

m3 <- dt[cerad ==0, clm(f4, data = .SD)]
summary(m3)

fm4_null <- formula(paste0(c('factor(braak56) ~ I(cerad >2)', covars), 
                           collapse = ' + '))
fm4 <- formula(paste0(c('factor(braak56) ~ snp*I(cerad >2)', covars),
                      collapse = ' + '))
m4_null <- dt[, clm(fm4_null, data = .SD)]
m4 <- dt[, clm(fm4, data = .SD)]
summary(m4)
anova(m4_null, m4)

fm5_null <- formula(paste0(c('cerad_ord ~ braak_ord', covars), collapse = ' + '))
fm5 <- formula(paste0(c('cerad_ord ~ snp*braak_ord', covars), collapse = ' + '))
m5_null <- dt[, clm(fm5_null, data = .SD)]
m5 <- dt[, clm(fm5, data = .SD)]
summary(m5)
anova(m5_null, m5)
