#!/bin/bash

# Find if input variants are eQTLs in GTEx
# and write summary information to file
pheno=$1
results_file=$2
snp_col=$3
pval_col=$4

# 1. Take plink .assoc.logistic-style input file, then create list of 
# SNP with P < 1e-05
qtls_file=tmp/qtls_"$pheno".tmp
Rscript --vanilla code/gwas/write_snp_list.R \
  -f $results_file \
  --snp_col $snp_col \
  --pval_col $pval_col \
  --max_pval 0.00001 \
  -o tmp/"$pheno"_snps_p1e-5.txt
  
./code/coloc/01a_reformat_snps.R -f tmp/"$pheno"_snps_p1e-5.txt | \
  python ./code/coloc/find_qtls.py -p 22 -o $qtls_file

# # 2. Create tidy QTL summary table
# if test -f "$qtls_file"; then
#     Rscript --vanilla code/coloc/01b_qtls.R \
#       tmp/qtls_"$pheno".tmp \
#       tmp/rsid_key_"$pheno".tmp \
#       $pheno
# fi
