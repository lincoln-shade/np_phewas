
phenos=$(cat data/adni/adni_np.pheno | awk 'NR==1{for (i=3; i<NF; i++) printf $i "\n"; print $NF}')

for pheno in $phenos
do
python code/analysis/plink_logistic.py \
  -b data/adni/adni_np \
  -o output/gwas/adni/"$pheno" \
  -p data/adni/adni_np.pheno \
  -c data/adni/adni_np.covar \
  --pheno_name "$pheno" &
done

for pheno in $phenos 
do 
echo tail output/gwas/adni/"$pheno".assoc.logistic
done