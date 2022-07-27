
phenos=$(cat data/mega/mega_np.pheno | awk 'NR==1{for (i=3; i<(NF); i++) printf $i "\n"; print $NF}')

# for pheno in $phenos
# do
# python code/gwas/plink2_logistic.py \
#   -f data/mega/mega_np \
#   -o output/gwas/mega/apoe_chr19/mega \
#   -p data/mega/mega_np.pheno \
#   -c data/mega/mega_np_apoe.covar \
#   -r data/mega/related_rm/"$pheno".remove \
#   --chrom 19 \
#   --pheno_name "$pheno"
# done

# for pheno in $phenos
# do
#   cat output/gwas/mega/apoe_chr19/mega."$pheno".glm.logistic | \
#     awk '$14 < 0.01 {print $3}' \
#     > tmp/"$pheno"_apoe_p0.01.txt
#   
#   plink2 \
#     --pfile data/mega/mega_np \
#     --extract tmp/"$pheno"_apoe_p0.01.txt \
#     --export A \
#     --export-allele data/mega/varID_ALT.txt \
#     --out tmp/"$pheno"_apoe_p0.01
# done

run_ord_reg () {
  ord_pheno=$1
  bin_pheno=$2
  ./code/gwas/ordinal/clm_regression.R \
  -c data/mega/mega_np_apoe.covar \
  -g tmp/"$bin_pheno"_apoe_p0.01.raw \
  -o tmp/"$ord_pheno"_apoe_ord_results.txt \
  -p data/mega/mega_np_ord.pheno \
  --phenotype "$ord_pheno" \
  -r data/mega/related_rm/"$ord_pheno".remove
}

# run_ord_reg caa_ord caa
# run_ord_reg cerad cerad3
# run_ord_reg braak braak56
run_ord_reg diffuse_abeta diffuse_abeta3
run_ord_reg late late23
run_ord_reg lewy_body lewy_body123
