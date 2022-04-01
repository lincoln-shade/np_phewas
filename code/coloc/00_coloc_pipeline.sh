#!/bin/bash
phenotype=$1

# the plink2 glm results for the phenotype
results=output/gwas/mega/mega."$phenotype".glm.logistic
# the prefix for the plink2 binary fileset used for analysis
data_prefix=data/mega/mega_np

bash code/coloc/01_qtls.sh \
  $phenotype \
  $results

top_qtls_file=tmp/top_qtls_"$phenotype".tmp

if test -f $top_qtls_file; then
    bash code/coloc/02_coloc.sh \
      $data_prefix \
      "$data_prefix".pheno \
      $results \
      output/coloc/mega/ \
      $phenotype
fi