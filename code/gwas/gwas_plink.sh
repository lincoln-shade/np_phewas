#!/bin/bash

phenotype=$1
out_file=output/gwas_results/adc_np_"$phenotype".assoc.logistic
echo CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P > $out_file

seq 1 22 | awk -v phe="$phenotype" '{print phe, $1}' | xargs -n 2 -P 22 bash code/analysis/gwas_plink_a.sh

for chr in $(seq 1 22)
  do
  cat output/tmp/adc_np_"$phenotype"_chr"$chr".assoc.logistic | awk 'NR>1{print $0}' >> $out_file
  done

rm output/tmp/adc_np_"$phenotype"_chr*