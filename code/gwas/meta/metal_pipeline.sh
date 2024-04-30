
meta_dir=./code/gwas/meta

# bash "$meta_dir"/prep_polmm_results_for_metal.sh
# bash "$meta_dir"/prep_saige_results_for_metal.sh

# for pheno in hs grossinf microinf microinf grossinf cerad braak caa arteriol atheroscler lewy_body late braak_apoe cerad_apoe late_apoe caa_apoe
#   do
#   metal "$meta_dir"/"$pheno"_metal.txt
#   done

# for pheno in braak cerad caa hs microinf grossinf late lewy_body arteriol atheroscler diffuse_abeta
#   do
#   Rscript --vanilla code/gwas/format_metal_results.R \
#     --metal_file output/gwas/metal/results/"$pheno"1.tbl \
#     --bim_file data/adc/adc_np.bim \
#     --output_file output/gwas/metal/results/"$pheno"1.csv 
#   done

for pheno in diffuse_abeta
do
python ./code/gwas/clump_snps.py \
  -b data/adc/adc_np \
  -r tmp/"$pheno"_magma_input.csv \
  -o output/gwas/metal/results/"$pheno"1_clump
done