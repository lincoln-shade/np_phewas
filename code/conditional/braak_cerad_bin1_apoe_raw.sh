
plink \
  --bfile data/mega/mega_np \
  --extract data/mega/coloc/apoe_snps.txt \
  --recode A \
  --out data/mega/conditional/apoe

plink \
  --bfile data/mega/mega_np \
  --extract data/mega/coloc/bin1_snps.txt \
  --recode A \
  --out data/mega/conditional/bin1