#!/bin/bash

# create exlusion files for each phenotype
cat data/plink/adc_np.pheno | awk 'NR==1{for (i=3; i<=NF; i++) print $i}' | xargs -n 1 Rscript --vanilla code/prep/11a_create_remove_file.R data/plink/adc_np.genome