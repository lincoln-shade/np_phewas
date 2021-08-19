#!/bin/bash

# Find if input variants are eQTLs in GTEx
# and write summary information to file

# 1. Take plink .assoc.logistic-style input file, then create list of SNP with P < 1e-05
cat $1 | awk '$12<0.00001 {print$2}' >> data/tmp/top_snps.tmp
# reformat variant names
Rscript code/coloc/01a_reformat_snps.R data/tmp/top_snps.tmp

# search for QTLs
grep --file=data/tmp/snps.tmp -r --include \*.signif_pairs.txt /data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_Eur/eqtls >> data/tmp/qtls.tmp
grep --file=data/tmp/snps.tmp -r --include \*.signif_pairs.txt /data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_sQTL_Eur/ >> data/tmp/qtls.tmp
# tidy QTL summary table
Rscript code/coloc/01b_qtls.R data/tmp/qtls.tmp data/tmp/rsid_key.tmp

# rm data/tmp/snps.tmp data/tmp/rsid_key.tmp data/tmp/qtls.tmp
