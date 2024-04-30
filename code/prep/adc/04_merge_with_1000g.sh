#!/bin/bash

#================================
# Step 3: Remove ethnic outliers
#================================
prefix=tmp/adc_qc1.tmp

# # convert to PLINK 1.9 fileset
# # (merging in PLINK 2.0 is currently limited in features)
# plink2 --pfile "$prefix" --make-bed --out "$prefix"
# 
# # merge adc and 1000g filesets
# python code/prep/merge_plink_filesets.py \
#   -b "$prefix" \
#   -b1 data/1000g/1000g \
#   -o tmp/adc_1000g_merged

# prune
plink \
  --bfile tmp/adc_1000g_merged \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/adc_1000g_merged_prune

# pca
plink \
  --bfile tmp/adc_1000g_merged \
  --no-pheno \
  --extract tmp/adc_1000g_merged_prune.prune.in \
  --pca 2 \
  --out tmp/adc_1000g_merged_pca
  