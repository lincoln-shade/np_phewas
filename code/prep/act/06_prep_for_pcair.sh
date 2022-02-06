#!/bin/bash

#=====================================
# prep for PC-AiR on related data set
#=====================================
prefix=act
prefix_out="$prefix"_np
# keep only european subjects
plink \
  --bfile tmp/"$prefix"_qc1.tmp \
  --keep tmp/"$prefix"_qced.txt \
  --make-bed \
  --out tmp/"$prefix"_qc3.tmp

# one final pass at variant QC
plink \
  --bfile tmp/"$prefix"_qc3.tmp \
  --maf 0.05 \
  --geno 0.05 \
  --hwe 1e-6 midp include-nonctrl \
  --make-bed \
  --out data/"$prefix"/"$prefix_out"

cat data/"$prefix"/"$prefix_out".fam | awk '{print 2, $2, 0, 0, 0, -9}' > data/"$prefix"/"$prefix_out".fam.tmp
mv data/"$prefix"/"$prefix_out".fam.tmp data/"$prefix"/"$prefix_out".fam

# create pruned dataset for PC-AiR
# prune
plink \
  --bfile data/"$prefix"/"$prefix_out" \
  --indep-pairwise 15000 1500 0.2 \
  --out data/"$prefix"/"$prefix_out"

plink \
  --bfile data/"$prefix"/"$prefix_out" \
  --extract data/"$prefix"/"$prefix_out".prune.in \
  --make-bed \
  --out data/"$prefix"/"$prefix_out"_pruned

  
# king kinship estimation
king -b data/"$prefix"/"$prefix_out".bed --kinship --prefix data/"$prefix"/"$prefix_out"
