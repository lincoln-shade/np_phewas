
library(data.table)
library(magrittr)
library(MASS)
library(mediation)
library(ordinal)

pheno <- setDT(readRDS('data/mega/mega_np.Rds'))[
  , .(FID, IID, late, braak)
]
covar <- fread('data/mega/mega_np.covar')
tmem_snps <- fread('data/mega/conditional/tmem106b.raw')
tmem_top_late_snp <- 'rs13237715'
tmem_top_snp_col <- colnames(tmem_snps)[
  grep(tmem_top_late_snp, colnames(tmem_snps))
]

covars <- colnames(covar)[3:ncol(covar)]

dt <- Reduce(
  merge, 
  list(pheno, 
       covar, 
       tmem_snps[, .(FID, IID, 'snp' = get(..tmem_top_snp_col))]
  )
)

dt_full <- copy(dt)
dt <- dt[complete.cases(dt)]

make_formula <- function(v, covars) {
  formula(paste(c(v, covars), collapse = " + "))
}
# mediation analysis
n_sims <- 1e4
f1 <- formula(paste0(c('braak ~ snp', covars), collapse = ' + '))
f2 <- formula(paste0(c('late ~ snp', covars), collapse = ' + '))
f3 <- formula(paste0(c('braak ~ snp', 'late', covars), collapse = ' + '))
iv_to_dv <- dt[, lm(f1, data = .SD)]
iv_to_m <- dt[, lm(f2, data = .SD)]
iv_m_to_dv <- dt[, lm(f3, data = .SD)]

mediation_linear <- mediate(iv_to_m, 
                          iv_m_to_dv, 
                          mediator = 'late', 
                          treat = 'snp', 
                          boot = FALSE, 
                          sims = n_sims)

saveRDS(mediation_linear, file = "output/mediation/mediation_linear.Rds")

# ordinal mediation analysis
dt[, late := ordered(late)]
dt[, braak := ordered(braak)]
iv_to_dv <- dt[, polr(f1, data = .SD, Hess = TRUE)]
iv_to_m <- dt[, polr(f2, data = .SD, Hess = TRUE)]
iv_m_to_dv <- dt[, polr(f3, data = .SD, Hess = TRUE)]

mediation_ordinal <- mediate(iv_to_m, 
                             iv_m_to_dv, 
                             mediator = 'late', 
                             treat = 'snp', 
                             boot = TRUE, 
                             sims = n_sims)

saveRDS(mediation_ordinal, file = "output/mediation/mediation_ordinal.Rds")

#---------------------------
# sensitivity analyses
#---------------------------
# redo basic analysis with clm function to make sure results are similar to 
# Stage 3 GWAS
dt_full[, braak_ord := as.ordered(braak)]
f1 <- make_formula("braak_ord ~ snp", covars)
m_braak <- dt_full[, clm(f1,
                         data = .SD)
] # results are similar, OR = 1.09, p = 0.00274

# remove participants without LATE-NC data
m_braak_no_lateNA <- dt_full[
  !is.na(late), 
  clm(f1, data = .SD)
] # OR = 1.17, p = 0.00419


# limit to Braak stage 1-4 to see if effect still there
m_braak1_4 <- dt_braak1_4[, clm(f1, data = .SD)]

dt_full[, braak := as.ordered(braak)]
m_braak <- dt_full[, clm(f1, data = .SD)]
m_braak1_4 <- dt_full[as.integer(braak) < 5, clm(f1, data = .SD)]

# adjust for LATE-NC
# OR = 1.27, p = 2.2e-5
m_braak_late <- dt_full[, clm(make_formula("braak_ord ~ snp + late", covars),
                              data = .SD)]

# limit to LATE-NC = 0 to see if effect still there
# OR = 1.22, 0.0063, N = 1283
m_late0 <- dt_full[as.integer(late) == 0, clm(f1, data = .SD)]

# perform analysis in individual cohorts
#        OR       P         N
# NACC   1.11     0.0045   5613
# ROSMAP 1.21     0.0010   1172
# ROS    1.27     0.026     576
# MAP    1.17     0.14      596
# ACT    0.95     0.63      614
# ADNI   4.27     0.033      39
m_braak_nacc <- dt_full[FID == "0", clm(f1, data = .SD)]
m_braak_rosmap <- dt_full[FID == "1", clm(f1, data = .SD)]
m_braak_ros <- dt_full[ROS == 1, clm(f1, data = .SD)]
m_braak_map <- dt_full[FID == "1" & ROS == 0, clm(f1, data = .SD)]
m_braak_act <- dt_full[FID == "2", clm(f1, data = .SD)]
m_braak_adni <- dt_full[FID == "ADNI", clm(f1, data = .SD)]
