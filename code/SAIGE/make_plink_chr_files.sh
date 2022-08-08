
# make chromosome-specific plink filesets
for chr in $(seq 1 22)
  do 
  plink \
    --bfile data/mega/mega_np \
    --chr "$chr" \
    --make-bed \
    --out data/mega/plink/mega_np_chr"$chr"
  done
