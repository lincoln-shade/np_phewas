
phenos=$(cat data/mega/mega_np_ord.pheno | cut -d " " -f 6- | awk 'NR == 1{print}')

for pheno in "vcid"
  do
  null_model=output/gwas/mega/polmm/"$pheno"_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/mega/mega_np.covar \
    --pheno data/mega/mega_np_ord.pheno \
    --plink tmp/mega_np_prune \
    --phenotype $pheno \
    --out $null_model \
    --ncores 100

  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/mega/SnpMatrix/SnpMatrix_list.Rds \
    -o output/gwas/mega/polmm/"$pheno"_polmm_results.Rds
  done



