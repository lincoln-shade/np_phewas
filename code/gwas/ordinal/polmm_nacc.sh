
mkdir output/gwas/adc/polmm

# create SNPMatrix
mkdir data/adc/SNPMatrix
./code/gwas/ordinal/polmm_create_SnpMatrix.R \
  --bim data/adc/adc_np.bim \
  --out data/adc/SNPMatrix/

phenos=$(cat data/mega/nacc_np.pheno | \
         cut -d " " -f 3- | \
         awk 'NR == 1{print}')

# main regression model
for pheno in diffuse_abeta arteriol atheroscler lewy_body caa cerad braak late 
  do
  null_model=output/gwas/adc/polmm/"$pheno"_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/adc/adc_np.covar \
    --pheno data/adc/adc_np_ord.txt \
    --plink data/adc/adc_np_pruned \
    --phenotype $pheno \
    --out $null_model \
    --ncores 100

  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/adc/SNPMatrix/SnpMatrix_list.Rds \
    -o output/gwas/adc/polmm/"$pheno"_polmm_results.Rds

  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r output/gwas/adc/polmm/"$pheno"_polmm_results.Rds \
  -o output/gwas/adc/polmm/"$pheno"_polmm_results.txt
  done

# run APOE-adjusted analysis for Chr19 SNPs
Rscript --vanilla ./code/gwas/ordinal/subset_snplist_chr19.R

for pheno in caa late braak cerad
  do
  null_model=output/gwas/adc/polmm/"$pheno"_apoe_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/adc/adc_np_apoe.covar \
    --pheno data/adc/adc_np_ord.txt \
    --plink data/adc/adc_np_pruned\
    --phenotype $pheno \
    --out $null_model \
    --ncores 100

  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/adc/SNPMatrix/SnpMatrix_list_chr19.Rds \
    -o output/gwas/adc/polmm/"$pheno"_apoe_polmm_results.Rds

  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r output/gwas/adc/polmm/"$pheno"_apoe_polmm_results.Rds \
  -o output/gwas/adc/polmm/"$pheno"_apoe_polmm_results.txt
  done

# cerad-adjusted (per reviewer request, not main manuscript)
for pheno in diffuse_abeta arteriol atheroscler lewy_body caa braak late
  do
  null_model=output/gwas/adc/polmm/"$pheno"_cerad_adj_null_model.Rds
  ./code/gwas/ordinal/polmm_null_model.R \
    --covar data/adc/adc_np_cerad.covar \
    --pheno data/adc/adc_np_ord.txt \
    --plink data/adc/adc_np_pruned \
    --phenotype $pheno \
    --out $null_model \
    --ncores 100

  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/adc/SNPMatrix/SnpMatrix_list.Rds \
    -o output/gwas/adc/polmm/"$pheno"_cerad_adj_polmm_results.Rds

  Rscript --vanilla ./code/gwas/rds_to_text.R \
  -r output/gwas/adc/polmm/"$pheno"_cerad_adj_polmm_results.Rds \
  -o output/gwas/adc/polmm/"$pheno"_cerad_adj_polmm_results.txt
  done