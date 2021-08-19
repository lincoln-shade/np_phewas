#==============================================================================
# Prep data for coloc analysis
#==============================================================================

library(pacman)
p_load(data.table, magrittr, stringi, rtracklayer, GEOquery, GENESIS, GWASTools, qusage, SeqArray, SNPRelate, arrow)

# cargs[1] = phenotype_id
# cargs[2] = tissue
# cargs[3] = chromosome
# cargs[4] = "eQTL" or "sQTL"
cargs <- commandArgs(trailingOnly = T)

phenotype_id_cargs <- cargs[1]
tissue_cargs <- cargs[2]
chr_cargs <- cargs[3]
qtl_type_cargs <- cargs[4]
# gwas data

load(paste0("data/tmp/chr", chr_cargs, "_gwas_sumstats.RData"))

#-------------------
# GTEx qtl data
#-------------------
if (qtl_type_cargs == "eQTL") {
  qtl <- as.data.table(read_parquet(paste0("/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/",
                      tissue_cargs,".v8.EUR.allpairs.chr", chr_cargs, ".parquet")))
}

if (qtl_type_cargs == "sQTL") {
  qtl <- as.data.table(read_parquet(paste0("/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_sQTL_all_associations/",
                      tissue_cargs,".v8.EUR.sqtl_allpairs.chr", chr_cargs, ".parquet")))
}


# phenotype ID
qtl <- qtl[phenotype_id == phenotype_id_cargs]

# split the variant_id column, make it into a data.table, and merge back
var.split <- stri_split_fixed(qtl[, variant_id], pattern = "_", n=5)
var.split <- as.data.table(t(matrix(do.call(c, var.split), nrow = 5)))
setnames(var.split, colnames(var.split), c("CHR", "BP_hg38", "A1", "A2", "build"))
var.split[, variant_id := paste(CHR, BP_hg38, A1, A2, build, sep = "_")]

qtl <- merge(qtl, var.split, "variant_id")
rm(var.split)

qtl[, BP_hg38 := as.integer(BP_hg38)]
qtl[, var.beta := slope_se**2]

tissues <- fread("/data_global/GTEx/gtex_tissue_sample_size.txt")
qtl[, n := tissues[tissue == tissue_cargs, round(n_genotyped*(1 - prop_non_eur))]] # number of GTEx participants who were in Artery_Aorta eQTL analysis
qtl[, type := "quant"] # QTL used linear regression
qtl <- qtl[, .(A1, BP_hg38, maf, pval_nominal, slope, var.beta, n, type)]
setnames(qtl, c("pval_nominal", "slope"), c("P", "beta"))

#------------------------------------
# Merge B-ASC and qtl summary stats
#------------------------------------

gwas_qtl <- merge(gwas, qtl, c("BP_hg38"), suffixes = c("_gwas", "_qtl"))

setorder(gwas_qtl, P_gwas)

fwrite(gwas_qtl, file= paste0("data/tmp/chr", chr_cargs, "_", phenotype_id_cargs, "_", tissue_cargs, "_gwas_qtl.txt"), sep=" ")

