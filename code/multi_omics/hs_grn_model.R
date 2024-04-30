library(data.table)
library(ordinal)
library(mediation)
library(MASS)
library(readxl)
rosmap_clin <- setDT(read_xlsx("/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))
rna <- fread("../np_phewas/data/rnaseq/rnaseq.txt.gz")
grn_col <- colnames(rna)[grep("ENSG00000030582", colnames(rna))]
rna_cols <- c("projid", "RIN", grn_col)
rna <- rna[, ..rna_cols]
rna[, projid := as.character(projid)]
dt_rna <- merge(rosmap_clin, rna, "projid")

# GRN and HS model
m_hs_grn = dt_rna[, glm(hspath_typ ~ ENSG00000030582.12 + msex + age_death + RIN + pmi, family = "binomial")]
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -11.34677    3.03992  -3.733 0.000190 ***
# ENSG00000030582.12   0.01027    0.02105   0.488 0.625673    
# msex                -0.10972    0.36263  -0.303 0.762215    
# age_death            0.09986    0.02770   3.605 0.000312 ***
# RIN                 -0.08210    0.16905  -0.486 0.627215    
# pmi                  0.01323    0.03167   0.418 0.676275    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# GRN and LATE model
dt_rna[, late := ordered(tdp_st4)]
m_late_grn = dt_rna[, clm(late ~ ENSG00000030582.12 + msex + age_death + RIN + pmi, data = .SD)]
summary(m_late_grn)
