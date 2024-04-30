# attempt to replicate results using Mayo CAA study data

library(data.table)
library(ordinal)

mccaa = fread("data/mc_caa/Metadata/MC-CAA_individual_human_metadata.csv")
mccaa_geno = fread("data/mc_caa/Metadata/MC-CAA_assay_snparray_metadata.csv")
dup_ids = c("18645", "18549", "815", "1287", "874", "3900", "925", "1014")
batchA_chr17 = fread("data/mc_caa/plink/mc_caa_batchA_chr17_rs72844606.raw", check.names = TRUE)
batchB_chr17 = fread("data/mc_caa/plink/mc_caa_batchB_chr17_rs72844606.raw", check.names = TRUE)
batchA_chr19 = fread("data/mc_caa/plink/mc_caa_batchA_chr19_rs7247551.raw", check.names = TRUE)
batchB_chr19 = fread("data/mc_caa/plink/mc_caa_batchB_chr19_rs7247551.raw", check.names = TRUE)

mccaa = mccaa[!(individualID %in% dup_ids)]
mccaa_geno[, individualID := as.character(specimenID)]
mc = merge(mccaa, mccaa_geno, "individualID")
mc[ageDeath == "90+", ageDeath := "90"] # ages >90 are truncated in synapse data
mc[, ageDeath := as.numeric(ageDeath)]
mc[, sex := factor(sex)]
mc[, apoeGenotype := factor(apoeGenotype)]
mc = mc[, .(individualID, AverageCAA, Braak, sex, ageDeath, PC1, PC2, PC3, PC4, PC5, apoeGenotype)]
mc[, individualID := as.integer(individualID)]

chr17 = rbind(batchA_chr17, batchB_chr17)
mc17 = merge(mc, chr17, by.x = "individualID", by.y = "IID")

model17 = mc17[Braak != 3.5, clm(as.ordered(Braak) ~ X17.8833591_T + ageDeath + sex + PC1 + PC2 + PC3 + PC4 + PC5)]

# major allele A is minor in Batch A; switch alleles
batchA_chr19[, X19.45454766_G := -X19.45454766_A + 2]
batchA_chr19[, X19.45454766_A := NULL]
chr19 = rbind(batchA_chr19, batchB_chr19)
mc19 = merge(mc, chr19, by.x = "individualID", by.y = "IID")

# Replication attempt model
model19 = mc19[, lm(sqrt(AverageCAA) ~ X19.45454766_G + apoeGenotype + ageDeath + sex + PC1 + PC2 + PC3)]
