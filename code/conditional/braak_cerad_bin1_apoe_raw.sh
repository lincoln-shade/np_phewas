
plink2 \
  --pfile data/mega/mega_np \
  --extract data/mega/coloc/apoe_snps.txt \
  --export A \
  --out data/mega/conditional/apoe

plink2 \
  --pfile data/mega/mega_np \
  --extract data/mega/coloc/bin1_snps.txt \
  --export A \
  --out data/mega/conditional/bin1