
study=act
# mkdir output/gwas/"$study"/polmm

# create SNPMatrix
# mkdir data/act/SNPMatrix
# ./code/gwas/ordinal/polmm_create_SnpMatrix.R \
#   --bim data/act/act_np.bim \
#   --out data/act/SNPMatrix/

# phenos=$(cat data/mega/"$study"_np_ord.pheno | \
#          cut -d " " -f 3- | \
#          awk 'NR == 1{print}')

for pheno in diffuse_abeta
  do
  # create null model
  null_model=output/gwas/"$study"/polmm/"$pheno"_null_model.Rds
  # ./code/gwas/ordinal/polmm_null_model.R \
  #   --covar data/"$study"/"$study"_np.covar \
  #   --pheno data/"$study"/"$study"_np.pheno \
  #   --plink data/"$study"/"$study"_np_pruned \
  #   --phenotype $pheno \
  #   --out $null_model \
  #   --ncores 100

  # run regression
  results=output/gwas/"$study"/polmm/"$pheno"_polmm_results.Rds
  ./code/gwas/ordinal/polmm_regression.R \
    -m $null_model \
    -l data/act/SNPMatrix/SnpMatrix_list.Rds \
    -o $results

  # Rscript --vanilla ./code/gwas/rds_to_text.R \
  # -r $results \
  # -o output/gwas/"$study"/polmm/"$pheno"_polmm_results.txt
  # 
  # python ./code/gwas/clump_snps.py \
  #   -r ./output/gwas/"$study"/polmm/"$pheno"_polmm_results.txt \
  #   -o output/gwas/"$study"/polmm/"$pheno"_polmm_results \
  #   -b data/act/act_np \
  #   --snp_field SNPID \
  #   --p_field pval.spa
  done

# for pheno in caa late braak cerad
#   do
#   null_model=output/gwas/act/polmm/"$pheno"_apoe_null_model.Rds
#   ./code/gwas/ordinal/polmm_null_model.R \
#     --covar data/act/act_np_apoe.covar \
#     --pheno data/act/act_np.pheno \
#     --plink data/act/act_np_pruned\
#     --phenotype $pheno \
#     --out $null_model \
#     --ncores 100
# 
#   ./code/gwas/ordinal/polmm_regression.R \
#     -m $null_model \
#     -l data/act/SNPMatrix/SnpMatrix_list_chr19.Rds \
#     -o output/gwas/act/polmm/"$pheno"_apoe_polmm_results.Rds
# 
#   Rscript --vanilla ./code/gwas/rds_to_text.R \
#   -r output/gwas/act/polmm/"$pheno"_apoe_polmm_results.Rds \
#   -o output/gwas/act/polmm/"$pheno"_apoe_polmm_results.txt
#   done

