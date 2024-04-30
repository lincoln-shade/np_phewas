library(data.table)
library(magrittr)
library(MASS)
library(mediation)
library(ordinal)

pheno <- setDT(readRDS('data/mega/mega_np.Rds'))[
  , .(FID, IID, caa_ord)
]

covar <- fread('data/mega/mega_np_apoe.covar')
caa_snp <- fread("data/mega/conditional/caa_apoe_snp.raw")
caa_snpid <- "rs4803779_C"
covars <- colnames(covar)[3:ncol(covar)]
dt <- Reduce(
  merge, 
  list(pheno, 
       covar, 
       caa_snp[, .(FID, IID, 'snp' = get(..caa_snpid))]
  )
)

