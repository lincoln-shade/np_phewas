#! /usr/bin/Rscript --vanilla
##-------------------------
## ordinal regression
##-------------------------

# use latest version of ordinal package (https://github.com/runehaubo/ordinal)
library(data.table)
library(magrittr)
library(ordinal)
library(stringi)
library(argparse)
library(dataPreparation)

parser <- ArgumentParser()
parser$add_argument("-c", "--covar", help="path to text file with covars")
parser$add_argument("-p", "--pheno", 
                    help="path to text file with ordinal phenotypes")
parser$add_argument("-g", "--geno", 
                    help="path to PLINK text file with genotypes")
parser$add_argument("-r", "--remove", 
                    help="path to text file with related ids to remove")
parser$add_argument("-o", "--out", 
                    help="path for output file")
parser$add_argument("--phenotype", help="phenotype name")
args <- parser$parse_args()

###############################################################################
# # test args
# args <- list(covar = "data/mega/mega_np_apoe.covar",
#              geno = "tmp/caa_apoe_p0.01.raw",
#              out = "tmp/caa_ord_apoe_ord_results.txt",
#              pheno = "data/mega/mega_np_ord.pheno",
#              phenotype = "caa_ord",
#              remove = "data/mega/related_rm/caa_ord.remove")
###############################################################################

raw <-fread(args$geno, header = T)
covar_data <- fread(args$covar, stringsAsFactors = TRUE)
phenotype_data <- fread(args$pheno, na.strings = "-1")
related_rm <- args$remove
phenotype <- args$phenotype
out <- args$out

# make sure FID and IID are character vectors for merging
covar_data[, `:=`(FID = as.character(FID), IID = as.character(IID))]
raw[, `:=`(FID = as.character(FID), IID = as.character(IID))]
phenotype_data[, `:=`(FID = as.character(FID), IID = as.character(IID))]
set(raw, 
    j = c("PAT", "MAT", "SEX", "PHENOTYPE"), 
    value = list(NULL, NULL, NULL, NULL)
)

# merge to create model matrix file
ord_data <- merge(phenotype_data[, c("FID", "IID", ..phenotype)], 
                  covar_data,
                  by = c("FID", "IID")
)
ord_data[, (phenotype) := as.ordered(get(..phenotype))]
ord_data <- ord_data[!is.na(get(phenotype))]

# remove individuals from related_rm folder file
try_fread <- function(path) {
  related_rm_ids <- try(fread(path, fill = T))
  if (inherits(related_rm_ids, "try-error")) {
    return(NULL)
  } else {
    return(related_rm_ids)
  }
}

related_rm_ids <- try_fread(related_rm)

if (!is.null(related_rm_ids)) {
  ord_data <- ord_data[!(IID %in% related_rm_ids$V2)]
}

# remove constant covariate columns
constant_col_indices <- which_are_constant(ord_data, keep_cols = 'FID')
constant_cols <- colnames(ord_data)[constant_col_indices]
if (length(constant_col_indices)) {ord_data[,  (constant_cols) := NULL]}
# merge with genotype data
ord_data <- merge(ord_data, raw, by = c("FID", "IID"))
rm(raw) # raw object very large, so remove after merging
gc()


# skip number of columns in covar_data + 1 
# (the phenotype) get to genetic variant columns
n_cols <- ncol(ord_data)
skip_cols <- ncol(covar_data) + 1L -length(constant_col_indices)
  
# replace colons with periods in SNP names, as colons in variable names 
# messes up regression
sub_colons <- function(x) { 
  x <- gsub(":", ".", x)
}

# add the character "v" in the front of all variant 
# names that do not start with "rs"
for (i in which(colnames(ord_data) %in% 
                colnames(ord_data)[-grep("rs", colnames(ord_data))
                                   ][-c(1:skip_cols)
                                     ])) {
  colnames(ord_data)[i] <- paste0("v", colnames(ord_data)[i])
}
  
colnames(ord_data) <- sub_colons(colnames(ord_data))



# initialize output table
results <- data.table(
  SNP = colnames(ord_data)[(skip_cols+1):n_cols],
  Beta = as.numeric(rep(NA, (n_cols - skip_cols))),
  SE = as.numeric(rep(NA, (n_cols - skip_cols))),
  P = as.numeric(rep(NA, (n_cols - skip_cols)))
)

# ordinal regression 
model_formula <- function(snp) {
  f <- formula(paste0(phenotype, 
                      " ~ ", 
                      snp,
                      " + ", 
                      paste(colnames(ord_data)[4:skip_cols],
                            collapse = " + ")))
}

run_clm <- function(i) {
  snp <- colnames(ord_data)[i]
  f <- model_formula(snp)
  m <- clm(f, data = ord_data)
  if (m$convergence$code == 0) {
    sm <- summary(m)
    set(
      results,
      i     = i - skip_cols,
      j     = c("Beta", "SE", "P"),
      value = list(
        sm$coefficients[snp, "Estimate"],
        sm$coefficients[snp, "Std. Error"],
        sm$coefficients[snp, "Pr(>|z|)"]
      )
    )
  }
}

for (i in (skip_cols + 1L):n_cols) {
  run_clm(i)
  if (i %% 100 == 0) {print(i)}
}
  
  
  
# round and format results
results[, `:=`(Beta = round(Beta, 4),
               SE = round(SE, 4),
               P = signif(P, 3))
]

# reformat rsIDs to match with those in pvar file
results[, A1 := stri_replace_first_regex(SNP, '.*_', '')]
results[, SNP := stri_replace_first_regex(SNP, '_[ACTG]*', '')]
results[, SNP := stri_replace_first_fixed(SNP, 'v', '')]
results[, SNP := stri_replace_all_fixed(SNP, '.', ':')]


# write outputs to file
fwrite(results, 
       file = out, 
       quote = F, 
       sep = " ",
       na = "NA"
)
