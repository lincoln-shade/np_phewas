library(data.table)
library(ordinal)
pheno <- fread("data/mega/mega_np_ord.pheno")
covar <- fread("data/mega/mega_np_age80.covar")
snp <- fread("tmp/abcc9_snps.raw")

covar[, FID := as.character(FID)]
dt <- merge(pheno[, .(FID, IID, arteriol)],
            covar,
            c("FID", "IID"))
dt <- merge(dt, snp[, .(FID, IID, rs1914361_G, rs704176_C)], c("FID", "IID"))

dt[, arteriol := as.ordered(arteriol)]
m1 <- dt[, clm(arteriol ~ rs1914361_G + age_death + msex + factor(FID))]
summary(m1)

m2 <- dt[, clm(arteriol ~ (rs1914361_G > 1) + age_death + msex + factor(FID))]
summary(m2)

m1 <- dt[, clm(arteriol ~ rs704176_C + age_death + msex + factor(FID))]
summary(m1)

m2 <- dt[, clm(arteriol ~ (rs704176_C > 1) + age_death + msex + factor(FID))]
summary(m2)
