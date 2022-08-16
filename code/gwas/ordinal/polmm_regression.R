#! /usr/bin/Rscript --vanilla
#=================================
# POLMM SNP regression
#=================================

library(data.table)
library(snpStats)
library(POLMM)
library(stringi)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-m", "--model", 
                    help = "path to POLMM null model .Rds file")
parser$add_argument("-l", "--snp_matrix_list",
                    default = "data/mega/SnpMatrix/SnpMatrix_list.Rds",
                    help = paste("path to .Rds file of vector", 
                                 "of SnpMatrix .Rds filenames"))
parser$add_argument("--snp_matrix", help = "path to SnpMatrix .Rds file")
parser$add_argument("-o", "--out", help = "path to .Rds file to save results")
args <- parser$parse_args()

if (is.null(args$model)) {
  stop("--model is a required argument")
}

if (!file.exists(args$model)) {stop(paste0(args$model, "doesn't exist"))}
if (!is.null(args$snp_matrix_list) & !file.exists(args$snp_matrix_list)) {
  stop(paste0(args$snp_matrix_list, "doesn't exist"))
}

if (!is.null(args$snp_matrix)) {
  snp_mat_files <- args$snp_matrix
  warning(paste0("Both --snp_matrix_list and --snp_matrix supplied.",
          " Using --snp_matrix and ignoring --snp_matrix_list."))
} else{
  snp_mat_files <- readRDS(args$snp_matrix_list)
}


null_model <- readRDS(args$model)
results_list <- vector("list", length = length(snp_mat_files))

for (i in 1:length(snp_mat_files)) {
  mat <- readRDS(snp_mat_files[i])
  g <- mat$genotypes@.Data
  class(g) <- "numeric"
  g <- g - 1
  g[g == -1] <- NA
  chrs <- mat$map$chromosome
  polmm <- setDT(POLMM(
    objNull = null_model,
    Geno.mtx = g,
    chrVec = chrs
    )
  )
  results_list[[i]] <- polmm
  rm(mat, g, chrs, polmm)
}

results <- rbindlist(results_list)

results[, `:=`(pval.spa = as.numeric(pval.spa),
               chr = as.integer(chr),
               missing.rate = as.numeric(missing.rate),
               Stat = as.numeric(Stat),
               VarW = as.numeric(VarW),
               VarP = as.numeric(VarP),
               pval.norm = as.numeric(pval.norm),
               MAF = as.numeric(MAF),
               beta = as.numeric(beta),
               switch.allele = as.logical(switch.allele))]

saveRDS(results, file = args$out)
