
phenos=$(cat data/mega/mega_np.pheno | awk 'NR==1{for (i=3; i<NF; i++) printf $i "\n"; print $NF}')

for pheno in part_pos1234 part_def1234
do
python code/gwas/plink2_logistic.py \
  -f data/mega/mega_np \
  -o output/gwas/mega/mega \
  -p data/mega/mega_np.pheno \
  -c data/mega/mega_np.covar \
  -r data/mega/related_rm/"$pheno".remove \
--pheno_name "$pheno"
done
