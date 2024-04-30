
# prep SAIGE results for METAL meta-analysis
nacc_res_dir=./output/gwas/adc/saige
rosmap_res_dir=./output/gwas/rosmap/saige
act_res_dir=./output/gwas/act/saige

metal_input_dir=./output/gwas/metal/input

for pheno in hs microinf grossinf
  do 
  
  nacc_results="$nacc_res_dir"/"$pheno"_saige_results.txt
  rosmap_results="$rosmap_res_dir"/"$pheno"_saige_results.txt
  act_results="$act_res_dir"/"$pheno"_saige_results.txt
  
  nacc_out="$metal_input_dir"/nacc_"$pheno".csv
  rosmap_out="$metal_input_dir"/rosmap_"$pheno".csv
  act_out="$metal_input_dir"/act_"$pheno".csv
  
  # Rscript --vanilla code/gwas/meta/prep_saige_results_for_metal.R \
  #   -r "$nacc_results" \
  #   --rm_alleles_rsid FALSE \
  #   -o "$nacc_out"
  #   
  # Rscript --vanilla code/gwas/meta/prep_saige_results_for_metal.R \
  #   -r "$rosmap_results" \
  #   --rm_alleles_rsid FALSE \
  #   -o "$rosmap_out"
    
    Rscript --vanilla code/gwas/meta/prep_saige_results_for_metal.R \
    -r "$act_results" \
    --rm_alleles_rsid TRUE \
    -o "$act_out"
  
  done
