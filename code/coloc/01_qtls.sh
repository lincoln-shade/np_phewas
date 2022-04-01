#!/bin/bash

# Find if input variants are eQTLs in GTEx
# and write summary information to file
pheno=$1
results_file=$2

# 1. Take plink .assoc.logistic-style input file, then create list of 
# SNP with P < 1e-05
qtls_file=tmp/qtls_"$pheno".tmp
cat $results_file | \
  awk '$14<0.00001 {print$3}' | \
  ./code/coloc/01a_reformat_snps.R $pheno | \
  python ./code/coloc/find_qtls.py -p 22 -o $qtls_file

# 2. Create tidy QTL summary table
if test -f "$qtls_file"; then
    Rscript --vanilla code/coloc/01b_qtls.R \
      tmp/qtls_"$pheno".tmp \
      tmp/rsid_key_"$pheno".tmp \
      $pheno
fi
