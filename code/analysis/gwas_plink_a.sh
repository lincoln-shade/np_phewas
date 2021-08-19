#!/bin/bash

phenotype=$1
chr=$2

plink \
  --bfile data/plink/adc_np \
  --logistic 'hide-covar' \
  --1 \
  --missing-phenotype -1 \
  --remove data/plink/remove/"$phenotype".remove \
  --pheno data/plink/adc_np.pheno \
  --pheno-name "$phenotype" \
  --covar data/plink/adc_np.covar \
  --ci 0.95 \
  --allow-no-sex \
  --chr "$chr" \
  --out output/tmp/adc_np_"$phenotype"_chr"$chr"
  