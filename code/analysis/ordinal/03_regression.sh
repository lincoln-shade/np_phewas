#!/bin/bash

phenotype=$1

<data/tmp/snp_list_index.tmp xargs -n 1 -P 22 \
  Rscript --vanilla --slave code/analysis/ordinal/03a_regression.R \
    data/plink/adc_np.covar data/adc_np_ord.txt data/related_rm/"$phenotype".remove $phenotype
