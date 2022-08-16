phenos=$(cat data/mega/mega_np.pheno | awk 'NR==1{for (i=3; i<NF; i++) printf $i "\n"; print $NF}')

for pheno in lewy_body123
do
awk 'NR==1{gsub("#CHROM", "CHR", $0) ; gsub("POS", "BP", $0) ; gsub("ID", "SNP", $0)}1' \
  output/gwas/mega/mega."$pheno".glm.logistic > tmp/mega."$pheno".p1.glm.logistic

./code/gwas/clump.sh \
  data/mega/mega_np \
  tmp/mega."$pheno".p1.glm.logistic \
  output/gwas/mega/mega."$pheno"
done

for pheno in lewy_body123
do
rm tmp/mega."$pheno".p1.glm.logistic
done