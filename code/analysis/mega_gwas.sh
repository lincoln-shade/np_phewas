
phenos=$(cat data/mega/mega_np.pheno | awk 'NR==1{for (i=3; i<NF; i++) printf $i "\n"; print $NF}')

# for pheno in $phenos
# do
# python code/analysis/plink_logistic.py \
#   -b data/mega/mega_np \
#   -o output/gwas/mega/"$pheno" \
#   -p data/mega/mega_np.pheno \
#   -c data/mega/mega_np.covar \
#   -r data/mega/related_rm/"$pheno".remove \
#   --pheno_name "$pheno" &
# done

for pheno in $phenos
do
./code/analysis/clump.sh data/mega/mega_np output/gwas/mega/"$pheno".assoc.logistic output/gwas/mega/"$pheno"_clumped
done