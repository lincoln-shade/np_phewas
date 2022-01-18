#!/bin/bash

#================================
# Step 3: Remove ethnic outliers
#================================
prefix=adc

# merge adc and 1000g filesets
python code/prep/merge_plink_filesets.py -b tmp/adc_qc1.tmp -b1 data/1000g/1000g -o tmp/adc_1000g_merged

# prune
plink \
  --bfile tmp/"$prefix"_1000g_merged \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/"$prefix"_1000g_merged_prune

# pca
plink \
  --bfile tmp/"$prefix"_1000g_merged \
  --no-pheno \
  --extract tmp/"$prefix"_1000g_merged_prune.prune.in \
  --pca 2 \
  --out tmp/"$prefix"_1000g_merged_pca
  