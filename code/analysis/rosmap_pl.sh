#!/bin/bash

phenotype=$1

plink \
  --bfile data/rosmap/rosmap_np \
  --logistic 'hide-covar' \
  --1 \
  --missing-phenotype -1 \
  --pheno data/rosmap/rosmap_np.pheno \
  --pheno-name "$phenotype" \
  --covar data/rosmap/rosmap_np.covar \
  --ci 0.95 \
  --allow-no-sex \
  --out output/gwas/rosmap/"$phenotype"
  