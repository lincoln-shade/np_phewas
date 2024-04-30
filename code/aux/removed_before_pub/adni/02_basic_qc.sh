#!/bin/bash

plink \
  --bfile /data_global/ADNI/ADNI_TOPMed_fromHohman_20211214/ADNI_NHW_imputed_final \
  --keep data/adni/adni_np_ids.txt \
  --make-bed \
  --out tmp/adni_qc0

cat tmp/adni_qc0.fam | awk '{print "ADNI", $2, $3, $4, $5, $6}' > tmp/adni_qc0_FID_ADNI.fam

# 20% variant missingness threshold
plink \
  --bfile tmp/adni_qc0 \
  --fam tmp/adni_qc0_FID_ADNI.fam \
  --geno 0.2 \
  --make-bed \
  --out tmp/adni_1a

# 20% IID missingness threshold
plink \
  --bfile tmp/adni_1a \
  --mind 0.2 \
  --make-bed \
  --out tmp/adni_1b

# 5% variant missingness threshold
plink \
  --bfile tmp/adni_1b \
  --geno 0.05 \
  --make-bed \
  --out tmp/adni_1c

# 5% IID missingness threshold
plink \
  --bfile tmp/adni_1c \
  --mind 0.05 \
  --make-bed \
  --out tmp/adni_1d

# create pruned SNP subset
# exclude regions of known high heterozygosity
plink \
  --bfile tmp/adni_1d \
  --exclude data/inversion.regions.txt \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/adni_1d_pruned

# create list of participants with unusually high heterozygosity
plink \
  --bfile tmp/adni_1d \
  --extract tmp/adni_1d_pruned.prune.in \
  --het \
  --out tmp/adni_1d

# check heterozygosity and create list of participants who fail QC
# (absolute heterozygosity rate > 3 sd from mean)
Rscript code/prep/check_heterozygosity.R -f tmp/adni_1d.het -o tmp/adni_fail_het.txt

# create PLINK fileset that keeps only those who pass heterozygosity QC
plink \
  --bfile tmp/adni_1d \
  --remove tmp/adni_fail_het.txt \
  --make-bed \
  --out tmp/adni_2
