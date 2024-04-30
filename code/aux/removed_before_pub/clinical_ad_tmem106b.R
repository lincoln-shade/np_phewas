#--------------------------------------------------------------------
# Investigate relationship between clinical Alzheimer's disease and 
# TMEM106B variants in study participants
#--------------------------------------------------------------------

# load neuropath, covariates, and TMEM106B variant data
library(data.table)
library(magrittr)
library(ordinal)
pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = "-1")
covar <- fread("data/mega/mega_np.covar")
tmem_snps <- fread('data/mega/conditional/tmem106b.raw')
tmem_top_late_snp <- 'rs13237715'
tmem_top_snp_col <- colnames(tmem_snps)[
  grep(tmem_top_late_snp, colnames(tmem_snps))
]

dt <- Reduce(
  merge, 
  list(pheno, 
       covar, 
       tmem_snps[, .(FID, IID, 'snp' = get(..tmem_top_snp_col))]
  )
)

dt <- dt[!is.na(braak)]

# rosmap clinical data
rosmap_cl <- fread("data/synapse/ROSMAP_clinical.csv")
rosmap_cl <- rosmap_cl[!is.na(cogdx), AD := ifelse(cogdx %in% 3:5, 1, 0)]
# nacc clinical data
na_strings <- c("-9", "-4", "9999", "999", "888")
nacc_uds <- fread("/data_global/nacc/investigator_nacc57.csv",
                  na.strings = na_strings) 
nacc_uds <- nacc_uds[
  NACCVNUM == NACCAVST][
    , .(NACCID, NACCALZD)][
      NACCALZD == 8, NACCALZD := 0
    ][
      , MDS := 0
    ]
setnames(nacc_uds, "NACCALZD", "AD")
# Minimal Data Set + NP
mds <- fread("/data_global/nacc/fardo09062019.csv", 
             na.strings = na_strings)
mds <- mds[, .(NACCID, CLINDEM)][CLINDEM == 2, CLINDEM := 0][, MDS := 1]
setnames(mds, "CLINDEM", "AD")

nacc <- rbind(nacc_uds, mds)
setorderv(nacc, cols = c("MDS"), order = -1L)
nacc <- nacc[!duplicated(NACCID)]

# merge rosmap and nacc
setnames(nacc, "NACCID", "IID")
setnames(rosmap_cl, "projid", "IID")
ad_pheno <- rbind(nacc[, .(IID, AD)], rosmap_cl[, .(IID, AD)])

# merge with data set
dt <- merge(dt, ad_pheno, by = "IID")

# analysis
make_formula <- function(v, covars) {
  formula(paste(c(v, covars), collapse = " + "))
}

# clinical AD diagnosis ~ TMEM106B
# OR = 0.98, P = 0.72
covars <- colnames(covar)[3:ncol(covar)]
f1 <- make_formula("AD ~ snp", covars)
m1 <- dt[, glm(f1, data = .SD, family = "binomial")]
m1.1 <- dt[!is.na(late), glm(f1, data = .SD, family = "binomial")]
summary(m1)

# Braak stage ~ TMEM106B
# OR = 1.11, P = 0.001
dt[, braak := ordered(braak)]
f2 <- make_formula("braak ~ snp", covars)
m2 <- dt[, clm(f2, data = .SD)]
summary(m2)

# LATE-NC ~ TMEM106B 
# OR = 0.71, P = 3.7e-8
dt[, late := ordered(late)]
f3 <- make_formula("late ~ snp", covars)
m3 <- dt[, clm(f3, data = .SD)]
summary(m3)
