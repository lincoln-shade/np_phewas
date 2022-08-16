#===================================================
# merge .pheno and .covar files for SAIGE analysis
#===================================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--pheno", "-p", help = "phenotype file")
parser$add_argument("--covar", "-c", help = "covariate file")
parser$add_argument("--out", "-o", help = "output file")
parser$add_argument("--na_string", default = "-1", help = "NA value string")
args <- parser$parse_args()

pheno <- fread(args$pheno, na.strings = args$na_string)
covar <- fread(args$covar)

pheno <- merge(pheno, covar)

fwrite(pheno, 
       file = args$out, 
       sep = " ",
       na = NA,
       quote = FALSE)
