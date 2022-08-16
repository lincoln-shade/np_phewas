cluster=vcid
pheno=PC1

covars=$(cat ./data/mega/mega_np.covar | \
         cut -d " " -f 3- | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

cat_covars=$(cat ./data/mega/mega_np.covar | \
         cut -d " " -f 3,5-21 | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

null_model_stage1=./output/multivar/"$cluster"_"$pheno"_stage1
Rscript --vanilla ./code/SAIGE/step1_fitNULLGLMM.R \
    --plinkFile=./tmp/mega_np_prune \
    --phenoFile=./data/mega/multivar/pcs_"$cluster".txt \
    --phenoCol="$pheno" \
    --covarColList="$covars" \
    --qCovarColList="$cat_covars" \
    --sampleIDColinphenoFile=IID \
    --nThreads=100 \
    --outputPrefix="$null_model_stage1" \
    --IsOverwriteVarianceRatioFile=TRUE \
    --traitType=quantitative

res_colname="$pheno"_res
res_pheno_file=./tmp/"$cluster"_"$res_colname".pheno
Rscript --vanilla ./code/SAIGE/saige_residuals.R \
  --model "$null_model_stage1".rda \
  --covar data/mega/mega_np.covar \
  --colname $res_colname \
  --sampleID_col IID \
  --output $res_pheno_file
  
null_model_stage2=./output/multivar/"$cluster"_"$pheno"_stage2
Rscript --vanilla ./code/SAIGE/step1_fitNULLGLMM.R \
    --plinkFile=./tmp/mega_np_prune \
    --phenoFile="$res_pheno_file" \
    --phenoCol="$res_colname" \
    --covarColList="$covars" \
    --qCovarColList="$cat_covars" \
    --sampleIDColinphenoFile=IID \
    --nThreads=100 \
    --outputPrefix="$null_model_stage2" \
    --IsOverwriteVarianceRatioFile=TRUE \
    --traitType=quantitative \
    --invNormalize=TRUE

saige_results=./output/multivar/"$cluster"_"$pheno"_saige_results_chr
for chr in $(seq 1 22)
    do
    plink=./data/mega/plink/mega_np_chr"$chr"
    Rscript --vanilla ./code/SAIGE/step2_SPAtests.R \
      --bedFile="$plink".bed \
      --bimFile="$plink".bim \
      --famFile="$plink".fam \
      --AlleleOrder=ref-first \
      --SAIGEOutputFile="$saige_results""$chr".txt \
      --chrom="$chr" \
      --GMMATmodelFile="$null_model_stage2".rda \
      --varianceRatioFile="$null_model_stage2".varianceRatio.txt \
      --is_output_moreDetails=TRUE \
      --LOCO=TRUE
    done

Rscript --vanilla ./code/SAIGE/bind_results.R \
  -p $saige_results