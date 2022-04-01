#!/bin/bash

# plink --bfile data/mega/mega_np_pruned --genome --min 0.15 --out data/mega/mega_np

# create exlusion files for each binary phenotype
cat data/mega/mega_np.pheno | \
  awk 'NR==1{for (i=3; i<=NF; i++) print $i}' | \
  xargs -n 1 \
  Rscript --vanilla code/prep/mega/04a_create_remove_file.R \
    data/mega/mega_np.genome \
    data/mega/mega_np.pheno \
    data/mega/related_rm