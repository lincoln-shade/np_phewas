# TMEM106B and GRN expression -> HS and TDP-43 vs other phenotypes

library(data.table)
library(readxl)
library(ordinal)
library(MASS)

np <- setDT(read_xlsx("/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))
rosmap <- fread("data/rosmap_data.csv")
rna_seq <- readRDS("data/rna_seq/rna_seq.Rds")
rna_seq[, projid := as.character(projid)]
tmem106b_col <- colnames(rna_seq)[grep("ENSG00000106460", colnames(rna_seq))]
grn_col <- colnames(rna_seq)[grep("ENSG00000030582", colnames(rna_seq))]
rna_seq_cols <- c("projid", "RIN", tmem106b_col, grn_col)
dt <- merge(np, rna_seq[, ..rna_seq_cols])
dt[, tdp_st4 := as.ordered(tdp_st4)]
dt[, ceradsc := 4 - ceradsc] # cerad scored as 1,2,3,4 severe -> none in ROSMAP
dt[, ceradsc := as.ordered(ceradsc)]
dt[, braaksc := as.ordered(braaksc)]
dt[, apoe := factor(apoe, 
                    levels = c("apoe33", "apoe34", "apoe44", 
                               "apoe23", "apoe22"))]
# TMEM106B models
# TDP-43 ~ TMEM106B
m1 <- dt[, clm(tdp_st4 ~ ENSG00000106460.13 + braaksc + msex + age_death + RIN + pmi)]

# HS ~ TMEM106B
m2 <- dt[, glm(hspath_typ ~ ENSG00000106460.13 + braaksc + msex + age_death + RIN + pmi,
               family = "binomial")]

# CERAD ~ TMEM106B

m3 <- dt[, clm(ceradsc ~ ENSG00000106460.13 + msex + age_death + RIN + pmi)]

# Braak ~ TMEM106B

m4 <- dt[, clm(braaksc ~ ENSG00000106460.13 + msex + age_death + RIN + pmi)]
m5 <- dt[, clm(braaksc ~ ENSG00000106460.13 + tdp_st4 + msex + age_death + RIN + pmi)]

# NFT burden ~ TMEM106B
m6 <- dt[, lm(sqrt(nft) ~ ENSG00000106460.13 + msex + age_death + RIN + pmi)]


# GRN models
grn_m1 <- dt[, clm(tdp_st4 ~ ENSG00000030582.12 + msex + age_death + RIN + pmi)]
grn_m2 <- dt[, clm(tdp_st4 ~ ENSG00000030582.12 + braaksc + msex + age_death + RIN + pmi)]
