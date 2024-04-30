
# find duplicates between NACC data set and Mayo CAA study
# 3 duplicates found with batch A
king \
  -b data/adc/adc_np.bed,data/mc_caa/Genetic_Variants/GWAS\ Genotypes\ -\ Batch\ A/MC-CAA_snpArray_batcha.bed \
  --duplicate \
  --prefix data/mc_caa/batchA_adc_dups

# 2 duplicates found in Batch B
king \
    -b data/adc/adc_np.bed,data/mc_caa/Genetic_Variants/GWAS\ Genotypes\ -\ Batch\ B/MC-CAA_snpArray_BatchB.bed \
    --duplicate \
    --prefix data/mc_caa/batchB_adc_dups

# 2 duplicates between batches A and B
king \
    -b data/mc_caa/Genetic_Variants/GWAS\ Genotypes\ -\ Batch\ A/MC-CAA_snpArray_batcha.bed,data/mc_caa/Genetic_Variants/GWAS\ Genotypes\ -\ Batch\ B/MC-CAA_snpArray_BatchB.bed \
    --duplicate \
    --prefix data/mc_caa/batchA_batchB_dups

