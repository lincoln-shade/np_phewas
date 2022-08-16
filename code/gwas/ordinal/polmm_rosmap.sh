
mkdir output/gwas/rosmap/polmm

phenos=$(cat data/mega/rosmap_np_ord.pheno | \
         cut -d " " -f 3- | \
         awk 'NR == 1{print}')

for pheno in $phenos
  do
  # output top clumped SNPs from NACC analysis to SnpMatrix
  snps=tmp/"$pheno"_nacc_snps.txt
  cat output/gwas/adc/polmm/"$pheno"_polmm_results.clumped | \
    awk 'NR>1{print $3}' \
    > $snps

  snpmat_prefix=data/mega/SnpMatrix/nacc_"$pheno"
  ./code/gwas/ordinal/polmm_create_SnpMatrix.R \
    --snp_list $snps \
    -o $snpmat_prefix

  # create null model
  null_model=output/gwas/rosmap/polmm/"$pheno"_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/mega/rosmap_np.covar \
    --pheno data/mega/rosmap_np_ord.pheno \
    --plink tmp/mega_np_prune \
    --phenotype $pheno \
    --out $null_model \
    --ncores 100

  # run regression
  results=output/gwas/rosmap/polmm/"$pheno"_polmm_results.Rds
  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    --snp_matrix "$snpmat_prefix"SnpMatrix.Rds \
    -o $results

  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r $results \
  -o output/gwas/rosmap/polmm/"$pheno"_polmm_results.txt
  done
