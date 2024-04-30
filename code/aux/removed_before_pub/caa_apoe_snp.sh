
snp=rs4803779
allele=C
export_allele=data/mega/conditional/caa_apoe_snp.txt
echo "$snp" "$allele" > "$export_allele"

plink2 \
  --pfile data/mega/mega_np \
  --export A \
  --export-allele "$export_allele" \
  --snp "$snp" \
  --out data/mega/conditional/caa_apoe_snp