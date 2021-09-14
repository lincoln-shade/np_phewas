#!/bin/bash

cat data/plink/adc_np.pheno | awk 'NR==1{for (i=3; i<=NF; i++) print $i}' | xargs -n 1 bash code/analysis/gwas_plink.sh