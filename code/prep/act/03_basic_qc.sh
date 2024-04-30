#!/bin/bash

#==============================================================
# basic QC
#==============================================================
prefix=act
#--------------------------------------------------------------
# make plink fileset of only those with autopsy data
#--------------------------------------------------------------

# apply initial MAF filter of 0.05%
plink --bfile tmp/"$prefix"_no_dup_rsid --keep data/"$prefix"_np_ids.txt --maf 0.005 --make-bed --out tmp/"$prefix"_np
# remove duplicated variants
cat tmp/"$prefix"_np.bim | awk '{print $2}' | sort | uniq -d > tmp/dup_vars.tmp
plink --bfile tmp/"$prefix"_np --exclude tmp/dup_vars.tmp --make-bed --out tmp/"$prefix"_np_uniq_varid

#-------------------------------------------------------------
# Step 1: basic variant and individual quality filters
# remove participants with missing genotyping rate of > 0.05
# remove variants with MAF < 0.05
# remove variants that fail HWE test
# remove variants that have missing genotype rate > 0.05
#--------------------------------------------------------------

# 20% variant missingness threshold
plink \
  --bfile tmp/"$prefix"_np_uniq_varid \
  --geno 0.2 \
  --make-bed \
  --out tmp/"$prefix"_qc1_geno.2.tmp

# 20% IID missingness threshold
plink \
  --bfile tmp/"$prefix"_qc1_geno.2.tmp \
  --mind 0.2 \
  --make-bed \
  --out tmp/"$prefix"_qc1_mind.2.tmp

# 5% variant missingness threshold
plink \
  --bfile tmp/"$prefix"_qc1_mind.2.tmp \
  --geno 0.05 \
  --make-bed \
  --out tmp/"$prefix"_qc1_geno.05.tmp

# 5% IID missingness threshold
plink \
  --bfile tmp/"$prefix"_qc1_geno.05.tmp \
  --mind 0.05 \
  --make-bed \
  --out tmp/"$prefix"_qc1_mind.05.tmp

# create pruned SNP subset
# exclude regions of known high heterozygosity
plink \
  --bfile tmp/"$prefix"_qc1_mind.05.tmp \
  --exclude raw_data/inversion.regions.txt \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/"$prefix"_qc1_pruned.tmp

# create list of participants with unusually high heterozygosity
plink \
  --bfile tmp/"$prefix"_qc1_mind.05.tmp \
  --extract tmp/"$prefix"_qc1_pruned.tmp.prune.in \
  --het \
  --out tmp/"$prefix"_qc1_het.tmp

# check heterozygosity and create list of participants who fail QC
# (absolute heterozygosity rate > 3 sd from mean)
Rscript --vanilla code/prep/act/03b_check_heterozygosity.R \
  tmp/"$prefix"_qc1_het.tmp.het \
  tmp/"$prefix"_qc1_fail_het.tmp

# create PLINK fileset that keeps only those who pass heterozygosity QC
plink \
  --bfile tmp/"$prefix"_qc1_mind.05.tmp \
  --remove ./tmp/"$prefix"_qc1_fail_het.tmp \
  --make-bed \
  --out tmp/"$prefix"_qc1.tmp

