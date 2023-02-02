
library(data.table)
library(ordinal)
mcbb <- fread("data/mcbb/MayoRNAseq_individual_metadata_031422.csv")
snp <- fread("data/mcbb/tmem106b.raw", check.names = TRUE)

dt <- merge(mcbb, snp, by.x = "individualID", by.y = "IID")

dt[, Braak_rd := as.ordered(floor(Braak))]
dt[ageDeath == "90_or_over", ageDeath := "90"]
dt[, ageDeath := as.numeric(ageDeath)]
m <- dt[race == "White", clm(as.ordered(Braak_rd) ~ X7.12269593_A + sex + ageDeath + factor(apoeGenotype))]
summary(m)
