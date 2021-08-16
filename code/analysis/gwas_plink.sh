#!/bin/bash

phenotype=ASC

echo CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P > output/gwas_results/adc_np_"$phenotype".assoc.logistic

for chr in 1..22
  do
  plink \
    --bfile data/plink/adc_np \
    --logistic 'hide-covar' \
    --1 \
    --missing-phenotype -1 \
    --remove data/plink/remove/ASC.remove \
    --pheno data/plink/adc_np.pheno \
    --pheno-name "$phenotype" \
    --covar data/plink/adc_np.covar \
    --ci 0.95 \
    --allow-no-sex \
    --chr "$chr"
    --out output/tmp/adc_np_"$phenotype"_chr"$chr"
  done

cat output/tmp/adc_np_"$phenotype"_chr"$chr".assoc.logistic | awk 'NR>1{print $0}' >> output/gwas_results/adc_np_"$phenotype".assoc.logistic
