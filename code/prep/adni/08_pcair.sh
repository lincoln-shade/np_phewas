#!/bin/bash

#=====================================
# prep for PC-AiR on related data set
#=====================================
prefix=adni
prefix_out="$prefix"_np

# one final pass at variant QC
plink \
  --bfile tmp/"$prefix"_3 \
  --maf 0.05 \
  --geno 0.05 \
  --hwe 1e-6 midp include-nonctrl \
  --make-bed \
  --keep data/"$prefix"/ids_qced_not_dups.txt \
  --out data/"$prefix"/"$prefix_out"

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

code/prep/pcair.R \
  -b data/adni/adni_np_pruned \
  -k data/adni/adni_np.kin \
  -o data/adni/mypcair.Rdata