
mkdir output/gwas/rosmap/polmm

# create SNPMatrix
mkdir data/rosmap/SNPMatrix
./code/gwas/ordinal/polmm_create_SnpMatrix.R \
  --bim data/rosmap/rosmap_np.bim \
  --out data/rosmap/SNPMatrix/


phenos=$(cat data/mega/rosmap_np_ord.pheno | \
         cut -d " " -f 3- | \
         awk 'NR == 1{print}')


for pheno in cerad diffuse_abeta
  do
  # create null model
  null_model=output/gwas/rosmap/polmm/"$pheno"_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/rosmap/rosmap_np.covar \
    --pheno data/rosmap/rosmap_np.pheno \
    --plink data/rosmap/rosmap_np_pruned \
    --phenotype $pheno \
    --out $null_model \
    --ncores 100

  # run regression
  results=output/gwas/rosmap/polmm/"$pheno"_polmm_results.Rds
  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/rosmap/SNPMatrix/SnpMatrix_list.Rds \
    -o $results

  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r $results \
  -o output/gwas/rosmap/polmm/"$pheno"_polmm_results.txt

  python ./code/gwas/clump_snps.py \
    -r ./output/gwas/rosmap/polmm/"$pheno"_polmm_results.txt \
    -o output/gwas/rosmap/polmm/"$pheno"_polmm_results \
    -b data/rosmap/rosmap_np \
    --snp_field SNPID \
    --p_field pval.spa
  done

# run APOE-adjusted analysis for Chr19 SNPs
Rscript --vanilla ./code/gwas/ordinal/subset_snplist_chr19.R

for pheno in caa late braak cerad
  do
  null_model=output/gwas/rosmap/polmm/"$pheno"_apoe_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/rosmap/rosmap_np_apoe.covar \
    --pheno data/rosmap/rosmap_np.pheno \
    --plink data/rosmap/rosmap_np_pruned\
    --phenotype $pheno \
    --out $null_model \
    --ncores 100

  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/rosmap/SNPMatrix/SnpMatrix_list_chr19.Rds \
    -o output/gwas/rosmap/polmm/"$pheno"_apoe_polmm_results.Rds

  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r output/gwas/rosmap/polmm/"$pheno"_apoe_polmm_results.Rds \
  -o output/gwas/rosmap/polmm/"$pheno"_apoe_polmm_results.txt
  done
