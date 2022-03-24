
library(data.table)
library(magrittr)
library(MASS)
library(mediation)

pheno <- setDT(readRDS('data/mega/mega_np.Rds'))[
  , .(FID, IID, hs, late, late23)
]
covar <- fread('data/mega/mega_np.covar')
related_rm <- fread('data/mega/related_rm/hs.remove')
tmem_snps <- fread('data/mega/conditional/tmem106b.raw')
tmem_top_late_snp <- 'rs1060700'
tmem_top_snp_col <- colnames(tmem_snps)[grep(tmem_top_late_snp, colnames(tmem_snps))]

covars <- colnames(covar)[3:ncol(covar)] %>% 
  .[!(. %in% c('act', 'ACT2', 'ACT3', 'adni'))]

dt <- Reduce(merge, 
             list(pheno, 
                  covar, 
                  tmem_snps[, .(FID, IID, 'snp' = get(..tmem_top_snp_col))]))

dt <- dt[!(IID %in% related_rm$V2)] %>% 
  .[complete.cases(.)]

# binary mediation analysis
f1 <- formula(paste0(c('hs ~ snp', covars), collapse = ' + '))
f2 <- formula(paste0(c('late23 ~ snp', covars), collapse = ' + '))
f3 <- formula(paste0(c('hs ~ snp', 'late23', covars), collapse = ' + '))
iv_to_dv <- dt[, glm(f1, data = .SD, family = binomial())]
iv_to_m <- dt[, glm(f2, data = .SD, family = binomial())]
iv_m_to_dv <- dt[, glm(f3, data = .SD, family = binomial())]

mr_binary <- mediate(iv_to_m, 
                     iv_m_to_dv, 
                     mediator = 'late23', 
                     treat = 'snp', 
                     boot = TRUE, 
                     sims = 100)

# ordinal late mediation analysis
dt[, late_ord := ordered(late)]
f4 <- formula(paste0(c('late_ord ~ snp', covars), collapse = ' + '))
f5 <- formula(paste0(c('hs ~ snp', 'late_ord', covars), collapse = ' + '))
iv_to_m <- dt[, polr(f4, data = .SD, Hess = TRUE)]
iv_m_to_dv <- dt[, glm(f5, data = .SD, family = binomial())]

mr_ord <- mediate(iv_to_m, 
                     iv_m_to_dv, 
                     mediator = 'late_ord', 
                     treat = 'snp', 
                     boot = TRUE, 
                     sims = 200)

# ordinal hs -> late mediation analysis
dt[, late_ord := ordered(late)]
f4 <- formula(paste0(c('hs ~ snp', covars), collapse = ' + '))
f5 <- formula(paste0(c('late23 ~ snp', 'hs', covars), collapse = ' + '))
iv_to_m <- dt[, glm(f4, data = .SD, family = binomial())]
iv_m_to_dv <- dt[, glm(f5, data = .SD, family = binomial())]

mr_ord <- mediate(iv_to_m, 
                  iv_m_to_dv, 
                  mediator = 'hs', 
                  treat = 'snp', 
                  boot = TRUE, 
                  sims = 200)
