
pheno=hs

for cohort in nacc_old nacc_new nacc_b2000_NA_adni rosmap act
do
python code/gwas/plink2_logistic.py \
  -f data/mega/mega_np \
  -o output/aux/hs_snps/"$cohort" \
  -p data/mega/mega_np.pheno \
  -c data/mega/mega_np.covar \
  -r data/mega/related_rm/"$pheno".remove \
  --pheno_name "$pheno" \
  --extract data/aux/hs_snps/hs_2014_top_snps.txt \
  --keep data/aux/hs_snps/"$cohort"_ids.txt
done