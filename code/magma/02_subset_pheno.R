#! /use/bin/Rscript --vanilla
#=============================================
# Reformat and subset .pheno files for MAGMA
#=============================================

library(data.table)
library(argparse)

parser <- ArgumentParser()

parser$add_argument('-p', '--pheno', help  = 'path to phenotype file')
parser$add_argument('-r', '--remove', 
                    help  = 'path to file with sample IDs to exclude')
parser$add_argument('--phenotype', 
                    help  = 'phenotype name')
parser$add_argument('-o', '--out', help  = 'path to output file')


args <- parser$parse_args()

# #----------------------------------------------------
# # test args
# #----------------------------------------------------
# args$pheno <- 'data/mega/mega_np_ord.pheno'
# args$remove <- 'data/mega/related_rm/braak56.remove'
# args$phenotype <- 'braak56'
# args$out <- 'data/mega/magma/braak56.pheno'
# #----------------------------------------------------

phenotype <- args$phenotype
pheno <- fread(args$pheno, colClasses = list(character = c('FID', 'IID')),
               na.strings = '-1')
pheno <- pheno[!is.na(get(phenotype))]
rel_rm <- fread(args$remove, colClasses = c('character', 'character'))
setnames(rel_rm, colnames(rel_rm), colnames(pheno)[1:2])


pheno <- pheno[!rel_rm, on = names(rel_rm)]
pheno <- pheno[, .(FID, IID, get(..phenotype))]
pheno[, V3 := V3 + 1]
setnames(pheno, 'V3', paste(phenotype))

fwrite(pheno, file = args$out, sep = ' ', quote = FALSE)
