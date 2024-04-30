bfile=data/mega/mega_np

plink \
  --bfile "$bfile" \
  --exclude raw_data/inversion.regions.txt \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/mega_np_prune

plink \
  --bfile "$bfile" \
  --extract tmp/mega_np_prune.prune.in \
  --make-bed \
  --out tmp/mega_np_prune