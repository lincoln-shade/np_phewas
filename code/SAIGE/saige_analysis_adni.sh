#==================================================================
# Run SAIGE GWAS for binary variables: hs, microinf, and grossinf
#==================================================================
study=adni
mkdir ./output/gwas/$study/saige/

pheno_file=tmp/"$study"_saige.pheno
Rscript --vanilla ./code/SAIGE/make_pheno_file.R \
  --pheno data/mega/"$study"_np.pheno \
  --covar data/mega/"$study"_np.covar \
  --out $pheno_file

covars=$(cat $pheno_file | \
         cut -d " " -f 6- | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

cat_covars=$(cat $pheno_file | \
         cut -d " " -f 7,18- | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

phenos=$(cat $pheno_file | \
         cut -d " " -f 3-5 | \
         awk 'NR==1 {print $0}')
         


for pheno in $phenos
  do
  Rscript --vanilla ./code/SAIGE/step1_fitNULLGLMM.R \
    --plinkFile=./tmp/mega_np_prune \
    --phenoFile="$pheno_file" \
    --phenoCol="$pheno" \
    --covarColList="$covars" \
    --qCovarColList="$cat_covars" \
    --sampleIDColinphenoFile=IID \
    --nThreads=100 \
    --outputPrefix=./output/gwas/"$study"/saige/"$pheno"_saige \
    --IsOverwriteVarianceRatioFile=TRUE

  for chr in $(seq 1 22)
    do
    plink=./data/mega/plink/mega_np_chr"$chr"
    Rscript --vanilla ./code/SAIGE/step2_SPAtests.R \
      --bedFile="$plink".bed \
      --bimFile="$plink".bim \
      --famFile="$plink".fam \
      --AlleleOrder=ref-first \
      --SAIGEOutputFile=./output/gwas/"$study"/saige/"$pheno"_saige_results_chr"$chr".txt \
      --chrom="$chr" \
      --GMMATmodelFile=./output/gwas/"$study"/saige/"$pheno"_saige.rda \
      --varianceRatioFile=./output/gwas/"$study"/saige/"$pheno"_saige.varianceRatio.txt \
      --is_Firth_beta=TRUE \
      --pCutoffforFirth=0.05 \
      --is_output_moreDetails=TRUE \
      --LOCO=TRUE
    done

  # aggregate results into one file for each outcome
  Rscript --vanilla ./code/SAIGE/bind_results.R \
    -p ./output/gwas/"$study"/saige/"$pheno"_saige_results_chr
  done