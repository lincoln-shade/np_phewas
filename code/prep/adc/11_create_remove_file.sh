#!/bin/bash

# create exlusion files for each binary phenotype
cat data/adc/adc_np.pheno | \
  awk 'NR==1{for (i=3; i<=NF; i++) print $i}' | \
  xargs -n 1 Rscript --vanilla code/prep/adc/11a_create_remove_file.R data/adc/adc_np.genome data/adc/np_qced.Rds

# create exlusion files for each binary phenotype
cat data/adc/adc_np_ord.txt | \
  awk 'NR==1{for (i=3; i<=NF; i++) print $i}' | \
  xargs -n 1 Rscript --vanilla code/prep/adc/11a_create_remove_file.R data/adc/adc_np.genome data/adc/np_qced.Rds