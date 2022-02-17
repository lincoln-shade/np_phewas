#==============================================================================
# Prep data for coloc analysis
#==============================================================================

library(data.table)
library(magrittr)
library(stringi)
library(rtracklayer)
library(GEOquery)
library(GENESIS)
library(GWASTools)
library(qusage)
library(SeqArray)
library(SNPRelate)
library(arrow)

args <- commandArgs(trailingOnly = T)

gwas_phenotype <- args[1]
phenotype_id_args <- args[2] # QTL phenotype
tissue_args <- args[3]
chr_args <- args[4]
qtl_type_args <- args[5]
# gwas data

load(paste0("data/tmp/chr", chr_args, "_gwas_sumstats_", gwas_phenotype, ".RData"))

#-------------------
# GTEx qtl data
#-------------------
if (qtl_type_args == "eQTL") {
  qtl <- as.data.table(read_parquet(paste0("/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/",
                      tissue_args,".v8.EUR.allpairs.chr", chr_args, ".parquet")))
}

if (qtl_type_args == "sQTL") {
  qtl <- as.data.table(read_parquet(paste0("/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_sQTL_all_associations/",
                      tissue_args,".v8.EUR.sqtl_allpairs.chr", chr_args, ".parquet")))
}


# phenotype ID
qtl <- qtl[phenotype_id == phenotype_id_args]

# split the variant_id column, make it into a data.table, and merge back
var_split <- stri_split_fixed(qtl[, variant_id], pattern = "_", n=5)
var_split <- as.data.table(t(matrix(do.call(c, var_split), nrow = 5)))
setnames(var_split, colnames(var_split), c("CHR", "BP_hg38", "A1", "A2", "build"))
var_split[, variant_id := paste(CHR, BP_hg38, A1, A2, build, sep = "_")]

qtl <- merge(qtl, var_split, "variant_id")
rm(var_split)

qtl[, BP_hg38 := as.integer(BP_hg38)]
qtl[, var.beta := slope_se**2]

tissues <- fread("/data_global/GTEx/gtex_tissue_sample_size.txt")
qtl[, n := tissues[tissue == tissue_args, round(n_genotyped*(1 - prop_non_eur))]] # number of GTEx participants who were in Artery_Aorta eQTL analysis
qtl[, type := "quant"] # QTL used linear regression
qtl <- qtl[, .(A1, BP_hg38, maf, pval_nominal, slope, var.beta, n, type)]
setnames(qtl, c("pval_nominal", "slope"), c("P", "beta"))

#------------------------------------
# Merge B-ASC and qtl summary stats
#------------------------------------

gwas_qtl <- merge(gwas, qtl, c("BP_hg38"), suffixes = c("_gwas", "_qtl"))

setorder(gwas_qtl, P_gwas)

fwrite(gwas_qtl, 
       file = paste0(
           "data/tmp/chr", chr_args, "_", phenotype_id_args, "_", 
           tissue_args, "_gwas_qtl", gwas_phenotype, ".tmp"
       ), 
       sep=" "
)

