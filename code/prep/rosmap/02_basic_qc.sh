#!/bin/bash

plink \
  --bfile /data_global/ROSMAP/rosmap_topmed_20210303/ROSMAP_NHW_imputed_final \
  --fam /data_global/ROSMAP/rosmap_topmed_20210303/ROSMAP_NHW_imputed_final_converted.fam \
  --keep data/rosmap/np_ids.txt \
  --make-bed \
  --out tmp/rosmap_qc0

cat tmp/rosmap_qc0.fam | awk '{print 1, $2, $3, $4, $5, $6}' > tmp/rosmap_qc0_fid1.fam

# 5% variant missingness threshold
plink \
  --bfile tmp/rosmap_qc0 \
  --fam tmp/rosmap_qc0_fid1.fam \
  --geno 0.05 \
  --make-bed \
  --out tmp/rosmap_1a

# 5% IID missingness threshold
plink \
  --bfile tmp/rosmap_1a \
  --mind 0.05 \
  --make-bed \
  --out tmp/rosmap_1b

# filter variants by HWE test
# include non-controls because plink phenotypes aren't accurate
plink \
  --bfile tmp/rosmap_1b \
  --hwe 1e-6 midp include-nonctrl \
  --make-bed \
  --out tmp/rosmap_1c

# create pruned SNP subset
# exclude regions of known high heterozygosity
plink \
  --bfile tmp/rosmap_1c \
  --exclude data/inversion.regions.txt \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/rosmap_1c_pruned

# create list of participants with unusually high heterozygosity
plink \
  --bfile tmp/rosmap_1c \
  --extract tmp/rosmap_1c_pruned.prune.in \
  --het \
  --out tmp/rosmap_1c

# check heterozygosity and create list of participants who fail QC
# (absolute heterozygosity rate > 3 sd from mean)
Rscript code/prep/rosmap/02a_check_heterozygosity.R

# create PLINK fileset that keeps only those who pass heterozygosity QC
plink \
  --bfile tmp/rosmap_1c \
  --remove ./tmp/rosmap_1d_fail_het.txt \
  --make-bed \
  --out tmp/rosmap_2
