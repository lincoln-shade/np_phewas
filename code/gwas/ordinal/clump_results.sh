phenos=$(cat data/mega/mega_np_ord.pheno | awk 'NR==1{for (i=3; i<NF; i++) printf $i "\n"; print $NF}')

for pheno in $phenos
do
./code/gwas/clump.sh \
  data/mega/mega_np \
  output/gwas/mega/"$pheno"_ord_results.txt \
  output/gwas/mega/"$pheno"_ord_results
done
