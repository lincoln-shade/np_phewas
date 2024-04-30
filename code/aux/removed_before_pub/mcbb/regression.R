
library(data.table)
library(ordinal)
mcbb <- fread("data/mcbb/MayoRNAseq_individual_metadata_031422.csv")
snp <- fread("data/mcbb/tmem106b.raw", check.names = TRUE)
rnaseq <- fread("data/mcbb/rnaseq_gene_counts_normalized_tidy.csv")
tmem106b_col <- colnames(rnaseq)[grep("ENSG00000106460", colnames(rnaseq))]
rnaseq_cols <- c("specimenID", tmem106b_col) 
dt <- merge(mcbb, snp, by.x = "individualID", by.y = "IID")
dt <- merge(dt, rnaseq[, ..rnaseq_cols], by.x = "individualID", by.y = "specimenID")
dt[, Braak_rd := as.ordered(floor(Braak))]
dt[ageDeath == "90_or_over", ageDeath := "90"]
dt[, ageDeath := as.numeric(ageDeath)]
m <- dt[race == "White", clm(as.ordered(Braak_rd) ~ X7.12269593_A + sex + ageDeath + factor(apoeGenotype))]
summary(m)

m1 <- dt[race == "White", clm(as.ordered(Braak_rd) ~ log(ENSG00000106460) + sex + ageDeath +pmi)]
