#==================================================================
# Run SAIGE GWAS for binary variables: hs, microinf, and grossinf
#==================================================================
study=rosmap
mkdir ./output/gwas/$study/saige/

pheno_file=tmp/"$study"_cerad_adj_saige.pheno
Rscript --vanilla ./code/SAIGE/make_pheno_file.R \
  --pheno data/"$study"/"$study"_np.pheno \
  --covar data/"$study"/"$study"_np_cerad.covar \
  --out $pheno_file

covars=$(cat $pheno_file | \
         cut -d " " -f 14-,3 | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

cat_covars=$(cat $pheno_file | \
         cut -d " " -f 3,15,16 | \
         awk 'BEGIN {FS=" "; OFS = ",";} NR==1 {$1=$1; print $0}')

phenos=$(cat $pheno_file | \
         cut -d " " -f 4-6 | \
         awk 'NR==1 {print $0}')

for pheno in $phenos
  do
  Rscript --vanilla ./code/SAIGE/step1_fitNULLGLMM.R \
    --plinkFile=./data/rosmap/rosmap_np_pruned \
    --phenoFile="$pheno_file" \
    --phenoCol="$pheno" \
    --covarColList="$covars" \
    --qCovarColList="$cat_covars" \
    --sampleIDColinphenoFile=IID \
    --nThreads=100 \
    --outputPrefix=./output/gwas/"$study"/saige/"$pheno"_cerad_adj_saige \
    --IsOverwriteVarianceRatioFile=TRUE

  for chr in $(seq 1 22)
    do
    plink=./data/"$study"/plink/"$study"_np_chr"$chr"
    Rscript --vanilla ./code/SAIGE/step2_SPAtests.R \
      --bedFile="$plink".bed \
      --bimFile="$plink".bim \
      --famFile="$plink".fam \
      --AlleleOrder=ref-first \
      --SAIGEOutputFile=./output/gwas/"$study"/saige/"$pheno"_cerad_adj_saige_results_chr"$chr".txt \
      --chrom="$chr" \
      --GMMATmodelFile=./output/gwas/"$study"/saige/"$pheno"_cerad_adj_saige.rda \
      --varianceRatioFile=./output/gwas/"$study"/saige/"$pheno"_cerad_adj_saige.varianceRatio.txt \
      --is_Firth_beta=TRUE \
      --pCutoffforFirth=0.05 \
      --is_output_moreDetails=TRUE \
      --LOCO=TRUE
    done

  # aggregate results into one file for each outcome
  Rscript --vanilla ./code/SAIGE/bind_results.R \
    -p ./output/gwas/"$study"/saige/"$pheno"_cerad_adj_saige_results_chr
  done