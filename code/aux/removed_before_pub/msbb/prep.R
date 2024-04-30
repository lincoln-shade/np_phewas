# check for association between cg09555818 and cg13119609 with APOE and APOC2 methylation in MSBB
library(data.table)
library(magrittr)
library(stringi)

make_tidy <- function(m) {
  gene_id <- as.vector(m[1, ])
  specimenID <- rownames(m)[-1]
  m <- m[-1, ] %>% 
    apply(2, as.numeric) %>% 
    as.data.table()
  m[, specimenID := specimenID]
  setcolorder(m, "specimenID")
  
  colnames(m)[2:ncol(m)] <- gene_id
  m
}

caa_cpg = c("specimenID", "cg09555818", "cg13119609")
cpg_t = fread("data/msbb/ROU_13169_B01_CUS_MethylEPIC.illumina.beta.txt")
cpg = t(cpg_t)
cpg_tidy = make_tidy(cpg)
cpg_tidy = cpg_tidy[, ..caa_cpg]
cpg_tidy[, Sample_Name := as.integer(stri_replace_first_fixed(specimenID, "#", ""))]
cpg_tidy[, specimenID := NULL]
cpg_meta = fread("data/msbb/MSBB_assay_MethylationArray_metadata.csv")
cpg_tidy = merge(cpg_tidy, cpg_meta[, .(Sample_Name, specimenID)], "Sample_Name")
biospec_meta = fread("data/msbb/MSBB_biospecimen_metadata.csv")
cpg_tidy = merge(cpg_tidy, biospec_meta, "specimenID")

rna_t = fread("data/msbb/rnaseq_count_matrix_norm.tsv")
rna_tidy_cols = c("specimenID", "ENSG00000234906", "ENSG00000130203") # APOC2 and APOE
rna = t(rna_t)
rna_tidy = make_tidy(rna)
rna_tidy = rna_tidy[, ..rna_tidy_cols]
rna_meta = fread("data/msbb/MSBB_assay_rnaSeq_metadata.csv")
rna_tidy = merge(rna_tidy, rna_meta, "specimenID")
rna_tidy = merge(rna_tidy, biospec_meta, "specimenID")

rna_cpg_tidy = merge(
  rna_tidy[, .(individualID, ENSG00000234906, ENSG00000130203, cellType, RIN)],
  cpg_tidy[, .(individualID, cg09555818, cg13119609)],
  "individualID"
)

rna_cpg_tidy[, cor(cg09555818, ENSG00000130203), cellType]
rna_cpg_tidy[, cor(cg09555818, ENSG00000234906), cellType]

rna_cpg_tidy[cellType == "NeuN-_sox10-", summary(lm(ENSG00000130203 ~ cg13119609 + RIN))]

