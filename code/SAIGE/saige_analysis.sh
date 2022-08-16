#==================================================================
# Run SAIGE GWAS for binary variables: hs, microinf, and grossinf
#==================================================================

covars=$(cat ./data/mega/mega_np_pheno_covar.txt | \
         cut -d " " -f 6- | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

cat_covars=$(cat ./data/mega/mega_np_pheno_covar.txt | \
         cut -d " " -f 6,8-24 | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

phenos=$(cat ./data/mega/mega_np_pheno_covar.txt | \
         cut -d " " -f 3-5 | \
         awk 'NR==1 {print $0}')

for pheno in $phenos
  do
  Rscript --vanilla ./code/SAIGE/step1_fitNULLGLMM.R \
    --plinkFile=./tmp/mega_np_prune \
    --phenoFile=./data/mega/mega_np_pheno_covar.txt \
    --phenoCol="$pheno" \
    --covarColList="$covars" \
    --qCovarColList="$cat_covars" \
    --sampleIDColinphenoFile=IID \
    --nThreads=100 \
    --outputPrefix=./output/gwas/mega/saige/"$pheno"_saige \
    --IsOverwriteVarianceRatioFile=TRUE

  for chr in $(seq 1 22)
    do
    plink=./data/mega/plink/mega_np_chr"$chr"
    Rscript --vanilla ./code/SAIGE/step2_SPAtests.R \
      --bedFile="$plink".bed \
      --bimFile="$plink".bim \
      --famFile="$plink".fam \
      --AlleleOrder=ref-first \
      --SAIGEOutputFile=./output/gwas/mega/saige/"$pheno"_saige_results_chr"$chr".txt \
      --chrom="$chr" \
      --GMMATmodelFile=./output/gwas/mega/saige/"$pheno"_saige.rda \
      --varianceRatioFile=./output/gwas/mega/saige/"$pheno"_saige.varianceRatio.txt \
      --is_Firth_beta=TRUE \
      --pCutoffforFirth=0.05 \
      --is_output_moreDetails=TRUE \
      --LOCO=TRUE
    done
  done

# aggregate results into one file for each outcome
Rscript --vanilla ./code/SAIGE/cat_results.R

for pheno in $phenos
  do
  python ./code/gwas/clump_snps.py \
    -r output/gwas/mega/saige/"$pheno"_saige_results.txt \
    -o output/gwas/mega/saige/"$pheno"_saige_results \
    -b data/mega/mega_np \
    --snp_field MarkerID \
    --p_field p.value
  done