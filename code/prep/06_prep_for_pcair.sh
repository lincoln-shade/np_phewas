#!/bin/bash

#=====================================
# prep for PC-AiR on related data set
#=====================================
prefix=adc
prefix_out="$prefix"_np
# keep only european subjects
plink \
  --bfile data/tmp/"$prefix"_qc1.tmp \
  --keep data/"$prefix"_qced.txt \
  --make-bed \
  --out data/tmp/"$prefix"_qc3.tmp

# one final pass at variant QC
plink \
  --bfile data/tmp/"$prefix"_qc3.tmp \
  --maf 0.05 \
  --geno 0.05 \
  --hwe 1e-6 midp include-nonctrl \
  --make-bed \
  --out data/plink/"$prefix_out"

# create pruned dataset for PC-AiR
# prune
plink \
  --bfile data/plink/"$prefix_out" \
  --indep-pairwise 15000 1500 0.2 \
  --out data/plink/"$prefix_out"

plink \
  --bfile data/plink/"$prefix_out" \
  --extract data/plink/"$prefix_out".prune.in \
  --make-bed \
  --out data/plink/"$prefix_out"_pruned
  
# king kinship estimation
king -b data/plink/"$prefix_out".bed --kinship --prefix data/"$prefix_out"
