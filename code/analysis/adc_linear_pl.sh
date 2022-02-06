
phenotype=$1

plink \
  --bfile data/adc/adc_np \
  --pheno data/adc/adc_np_cont.pheno \
  --missing-phenotype -1 \
  --covar data/adc/adc_np.covar \
  --linear 'hide-covar' \
  --allow-no-sex \
  --ci 0.95 \
  --maf 0.05 \
  --pheno-name "$phenotype" \
  --remove data/adc/related_rm/"$phenotype".remove \
  --out output/gwas/adc/"$phenotype"