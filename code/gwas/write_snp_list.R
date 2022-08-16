#! /usr/bin/Rscript --vanilla
#==============================================
# Input text file with SNP IDs as one column
# and output list of SNP IDs
#==============================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-f", "--file", help = "text file SNPs as a column")
parser$add_argument("--snp_col", default = "SNP", help = "SNP column name")
parser$add_argument("-o", "--out", default = "./snp_list.txt",
                    help = "output filename")
args <- parser$parse_args()

snps <- fread(args$file)
fwrite(snps[, .(get(args$snp_col))],
       file = args$out,
       col.names = FALSE,
       quote = FALSE)

