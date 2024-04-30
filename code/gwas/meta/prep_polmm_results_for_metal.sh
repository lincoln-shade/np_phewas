
# prep POLMM results for METAL meta-analysis
nacc_res_dir=./output/gwas/adc/polmm
rosmap_res_dir=./output/gwas/rosmap/polmm
act_res_dir=./output/gwas/act/polmm

nacc_bim=./data/adc/adc_np.bim
rosmap_bim=./data/rosmap/rosmap_np.bim
act_bim=./data/act/act_np.bim

metal_input_dir=./output/gwas/metal/input

for pheno in diffuse_abeta
  do

  nacc_results="$nacc_res_dir"/"$pheno"_polmm_results.txt
  rosmap_results="$rosmap_res_dir"/"$pheno"_polmm_results.txt
  act_results="$act_res_dir"/"$pheno"_polmm_results.txt

  nacc_null="$nacc_res_dir"/"$pheno"_null_model.Rds
  rosmap_null="$rosmap_res_dir"/"$pheno"_null_model.Rds
  act_null="$act_res_dir"/"$pheno"_null_model.Rds

  nacc_out="$metal_input_dir"/nacc_"$pheno".csv
  rosmap_out="$metal_input_dir"/rosmap_"$pheno".csv
  act_out="$metal_input_dir"/act_"$pheno".csv

  Rscript --vanilla code/gwas/meta/prep_polmm_results_for_metal.R \
    -r "$nacc_results" \
    --bim "$nacc_bim" \
    -n "$nacc_null" \
    --rm_alleles_rsid FALSE \
    -o "$nacc_out"

  Rscript --vanilla code/gwas/meta/prep_polmm_results_for_metal.R \
    -r "$rosmap_results" \
    --bim "$rosmap_bim" \
    -n "$rosmap_null" \
    --rm_alleles_rsid FALSE \
    -o "$rosmap_out"
#     
#     Rscript --vanilla code/gwas/meta/prep_polmm_results_for_metal.R \
#     -r "$act_results" \
#     --bim "$act_bim" \
#     -n "$act_null" \
#     --rm_alleles_rsid TRUE \
#     -o "$act_out"
#   
  done

# apoe-adjusted
for pheno in braak cerad late
  do

  nacc_results="$nacc_res_dir"/"$pheno"_apoe_polmm_results.txt
  rosmap_results="$rosmap_res_dir"/"$pheno"_apoe_polmm_results.txt
  act_results="$act_res_dir"/"$pheno"_apoe_polmm_results.txt

  nacc_null="$nacc_res_dir"/"$pheno"_apoe_null_model.Rds
  rosmap_null="$rosmap_res_dir"/"$pheno"_apoe_null_model.Rds
  act_null="$act_res_dir"/"$pheno"_apoe_null_model.Rds

  nacc_out="$metal_input_dir"/nacc_"$pheno"_apoe.csv
  rosmap_out="$metal_input_dir"/rosmap_"$pheno"_apoe.csv
  act_out="$metal_input_dir"/act_"$pheno"_apoe.csv

  Rscript --vanilla code/gwas/meta/prep_polmm_results_for_metal.R \
    -r "$nacc_results" \
    --bim "$nacc_bim" \
    -n "$nacc_null" \
    --rm_alleles_rsid FALSE \
    -o "$nacc_out"

  Rscript --vanilla code/gwas/meta/prep_polmm_results_for_metal.R \
    -r "$rosmap_results" \
    --bim "$rosmap_bim" \
    -n "$rosmap_null" \
    --rm_alleles_rsid FALSE \
    -o "$rosmap_out"

    Rscript --vanilla code/gwas/meta/prep_polmm_results_for_metal.R \
    -r "$act_results" \
    --bim "$act_bim" \
    -n "$act_null" \
    --rm_alleles_rsid TRUE \
    -o "$act_out"

  done
