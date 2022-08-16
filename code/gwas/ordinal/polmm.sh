
phenos=$(cat data/mega/mega_np_ord.pheno | \
         cut -d " " -f 3- | \
         awk 'NR == 1{print}')

for pheno in $phenos
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

  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r output/gwas/mega/polmm/"$pheno"_polmm_results.Rds \
  -o output/gwas/mega/polmm/"$pheno"_polmm_results.txt

  python ./code/gwas/clump_snps.py \
    -r ./output/gwas/mega/polmm/"$pheno"_polmm_results.txt \
    -o output/gwas/mega/polmm/"$pheno"_polmm_results \
    -b data/mega/mega_np \
    --snp_field SNPID \
    --p_field pval.spa
  done

# run APOE-adjusted analysis for Chr19 SNPs
Rscript --vanilla ./code/gwas/ordinal/subset_snplist_chr19.R

for pheno in caa_ord
  do
  null_model=output/gwas/mega/polmm/"$pheno"_apoe_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/mega/mega_np_apoe.covar \
    --pheno data/mega/mega_np_ord.pheno \
    --plink tmp/mega_np_prune \
    --phenotype $pheno \
    --out $null_model \
    --ncores 100
  
  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/mega/SnpMatrix/SnpMatrix_list95-97.Rds \
    -o output/gwas/mega/polmm/"$pheno"_apoe_polmm_results.Rds
  
  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r output/gwas/mega/polmm/"$pheno"_apoe_polmm_results.Rds \
  -o output/gwas/mega/polmm/"$pheno"_apoe_polmm_results.txt
  
  python ./code/gwas/clump_snps.py \
    -r ./output/gwas/mega/polmm/"$pheno"_apoe_polmm_results.txt \
    -o output/gwas/mega/polmm/"$pheno"_apoe_polmm_results \
    -b data/mega/mega_np \
    --snp_field SNPID \
    --p_field pval.spa
  done