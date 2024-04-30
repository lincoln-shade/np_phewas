
chrom=7
bp=12236192
radius=200000

cat output/gwas/mega/mega.late23.glm.logistic | \
  awk -v CHR="$chrom" -v BP="$bp" -v RAD="$radius" '($1 == CHR) && ($2 >= BP-RAD) && ($2 <= BP+RAD){print $3, $6}' \
  > data/mega/conditional/tmem106b_snps.txt

plink2 \
  --pfile data/mega/mega_np \
  --export A \
  --export-allele data/mega/conditional/tmem106b_snps.txt \
  --extract data/mega/conditional/tmem106b_snps.txt \
  --out data/mega/conditional/tmem106b