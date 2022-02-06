
pref=data/act
for variable in $(cat data/act/act_pheno.txt)
do 
python code/analysis/plink_logistic.py -b "$pref"/act_np \
  -o output/gwas/act/"$variable" -r "$pref"/act_related.remove \
  -p "$pref"/act_np.pheno -c "$pref"/act_np.covar --pheno_name "$variable" &
done