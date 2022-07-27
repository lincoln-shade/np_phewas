#! /usr/bin/Rscript --vanilla
#======================================
# Fit null model for POLMM regression
#======================================

library(data.table)
library(POLMM)
library(argparse)
library(dataPreparation)
library(parallel)


parser <- ArgumentParser()
parser$add_argument("-c", "--covar", help="path to text file with covars")
parser$add_argument("-p", "--pheno", 
                    help="path to text file with ordinal phenotypes")
parser$add_argument("-k", "--grm", 
                    help="path to SparseGRM .Rds file. Not currently used.")
parser$add_argument("-b", "--plink", help="plink fileset prefix")
parser$add_argument("-o", "--out", 
                    help="path for output file")
parser$add_argument("--phenotype", help="phenotype name")
parser$add_argument("--ncores",
                    help = "number of CPU cores to use (default is 1/2 max)")
args <- parser$parse_args()

# #---------------------------------------------
# # test args
# args <- list(
#   covar = "data/mega/mega_np.covar",
#   grm = NULL,
#   out = "output/gwas/mega/polmm/diffuse_abeta_null_model.Rds",
#   pheno = "data/mega/mega_np_ord.pheno",
#   phenotype = "diffuse_abeta",
#   plink = "tmp/mega_np_prune",
#   ncores = "80"
# )
# #---------------------------------------------

covar_data <- fread(args$covar, stringsAsFactors = TRUE)
phenotype_data <- fread(args$pheno, na.strings = "-1")
phenotype <- args$phenotype
if (length(args$grm) > 0) {grm <- readRDS(args$grm)}
if (!is.null(args$ncores)) {
  args$ncores <- as.integer(args$ncores)
} else {
  args$ncores <- detectCores(all.tests = FALSE, logical = TRUE) / 2
  }

covar_data[, `:=`(FID = as.character(FID), IID = as.character(IID))]
phenotype_data[, `:=`(FID = as.character(FID), IID = as.character(IID))]
ord_data <- merge(phenotype_data[, c("FID", "IID", ..phenotype)], 
                  covar_data,
                  by = c("FID", "IID"))
ord_data[, (phenotype) := as.ordered(get(..phenotype))]
ord_data <- ord_data[!is.na(get(phenotype))]

# remove constant covariate columns
constant_col_indices <- which_are_constant(ord_data, keep_cols = 'FID')
constant_cols <- colnames(ord_data)[constant_col_indices]
if (length(constant_col_indices)) {ord_data[,  (constant_cols) := NULL]}

# create formula for null model
f <- formula(paste0(
  phenotype, " ~ ", paste(colnames(ord_data)[4:ncol(ord_data)], 
                          collapse = " + ")))

polmm_null_model <- 
  POLMM_Null_Model(
    formula = f,
    data = ord_data,
    subjData = ord_data$IID,
    PlinkFile = args$plink,
    control = list(numThreads = args$ncores)
  )

saveRDS(polmm_null_model, file = args$out)
