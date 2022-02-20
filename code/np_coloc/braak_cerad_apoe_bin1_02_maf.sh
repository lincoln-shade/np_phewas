

plink \
  --bfile data/mega/mega_np \
  --extract data/mega/coloc/apoe_snps.txt \
  --pheno data/mega/mega_np.pheno \
  --pheno-name braak56 \
  --missing-phenotype -1 \
  --1 \
  --remove data/mega/related_rm/braak56.remove \
  --freq \
  --allow-no-sex \
  --out data/mega/coloc/apoe

plink \
  --bfile data/mega/mega_np \
  --extract data/mega/coloc/bin1_snps.txt \
  --pheno data/mega/mega_np.pheno \
  --pheno-name braak56 \
  --missing-phenotype -1 \
  --1 \
  --remove data/mega/related_rm/braak56.remove \
  --freq \
  --allow-no-sex \
  --out data/mega/coloc/bin1