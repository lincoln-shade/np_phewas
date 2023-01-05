#=============================================================================
# Compare gene expression of COL4A1 and PIK3R5 with atherosclerosis and Braak
#=============================================================================

library(data.table)
library(ordinal)
library(readxl)
library(ggplot2)
rnaseq <- fread("data/rnaseq/rnaseq.txt.gz")
rnaseq[, projid := as.character(projid)]
# pheno <- fread("data/mega/mega_np_ord.pheno")
pheno <- read_xlsx("/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx")
setDT(pheno)
# covar <- fread("data/mega/mega_np.covar")
# dt <- merge(pheno, covar)
dt <- pheno[, .(projid, braaksc, cvda_4gp2, pmi, study, age_death, msex)]
col4a1 <- "ENSG00000187498"
pik3r5 <- "ENSG00000141506"
col4a1_col <- colnames(rnaseq)[grep(col4a1, colnames(rnaseq))]
pik3r5_col <- colnames(rnaseq)[grep(pik3r5, colnames(rnaseq))]
dt <- merge(dt, rnaseq[, c("projid", "RIN", ..col4a1_col, ..pik3r5_col)], by = "projid")

col4a1_lm_null <- dt[!is.na(cvda_4gp2), lm(log(ENSG00000187498.9) ~ msex + age_death + RIN +pmi,
                     data = .SD)]
col4a1_lm <- dt[, lm(log(ENSG00000187498.9) ~ ordered(cvda_4gp2) + msex + age_death + RIN +pmi,
                     data = .SD)]
anova(col4a1_lm_null, col4a1_lm)

pik3r5_lm_null <- dt[!is.na(braaksc), 
                     lm(log(ENSG00000141506.9) ~ msex + age_death + RIN +pmi,
                        data = .SD)]
pik3r5_lm <- dt[, lm(log(ENSG00000141506.9) ~ (braaksc > 3) + msex + age_death + RIN +pmi,
                     data = .SD)]
anova(pik3r5_lm_null, pik3r5_lm)

col4a1_clm <- dt[, clm(ordered(cvda_4gp2) ~ log(ENSG00000187498.9)+ study + msex + age_death + RIN +pmi,
                       data = .SD)]
pik3r5_clm <- dt[, clm(ordered(braaksc) ~ log(ENSG00000141506.9)+ study + msex + age_death + RIN +pmi,
                       data = .SD)]

dt[, ggplot(data=.SD, aes(braaksc > 3, log(ENSG00000141506.9)))] +
  geom_boxplot(aes(fill = braaksc > 3)) +
  geom_jitter(width = 0.2) +
  ylab("log(PIK3R5) expression") +
  xlab("Braak stage 4, 5, or 6") +
  theme_minimal() +
  theme(legend.position = "none")

  