#!/bin/bash

# Find if input variants are eQTLs in GTEx
# and write summary information to file
pheno=$1
results_file=output/gwas_results/adc_np_"$pheno".assoc.logistic

# 1. Take plink .assoc.logistic-style input file, then create list of SNP with P < 1e-05
qtls_file=data/tmp/qtls_"$pheno".tmp
cat $results_file | \
  awk '$12<0.00001 {print$2}' | \
  ./code/coloc/01a_reformat_snps.R $pheno | \
  python ./code/coloc/find_qtls.py -p 22 -o $qtls_file

# tidy QTL summary table
if test -f "$qtls_file"; then
    Rscript --vanilla code/coloc/01b_qtls.R data/tmp/qtls_"$pheno".tmp data/tmp/rsid_key_"$pheno".tmp $pheno
fi
# 
# rm data/tmp/qtls.tmp
