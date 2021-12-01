#!/bin/bash

#==============================================================
# basic QC
#==============================================================
prefix=act
#--------------------------------------------------------------
# make plink fileset of only those with autopsy data
#--------------------------------------------------------------

# apply initial MAF filter of 1%
plink --bfile data/tmp/"$prefix"_no_dup_rsid --keep data/"$prefix"_np_ids.txt --maf 0.01 --make-bed --out data/tmp/"$prefix"_np
# remove duplicated variants
cat data/tmp/"$prefix"_np.bim | awk '{print $2}' | sort | uniq -d > data/tmp/dup_vars.tmp
plink --bfile data/tmp/"$prefix"_np --exclude data/tmp/dup_vars.tmp --make-bed --out data/tmp/"$prefix"_np_uniq_varid

#-------------------------------------------------------------
# Step 1: basic variant and individual quality filters
# remove participants with missing genotyping rate of > 0.05
# remove variants with MAF < 0.05
# remove variants that fail HWE test
# remove variants that have missing genotype rate > 0.05
#--------------------------------------------------------------

# 20% variant missingness threshold
plink \
  --bfile data/tmp/"$prefix"_np_uniq_varid \
  --geno 0.2 \
  --make-bed \
  --out data/tmp/"$prefix"_qc1_geno.2.tmp

# 20% IID missingness threshold
plink \
  --bfile data/tmp/"$prefix"_qc1_geno.2.tmp \
  --mind 0.2 \
  --make-bed \
  --out data/tmp/"$prefix"_qc1_mind.2.tmp

# 5% variant missingness threshold
plink \
  --bfile data/tmp/"$prefix"_qc1_mind.2.tmp \
  --geno 0.05 \
  --make-bed \
  --out data/tmp/"$prefix"_qc1_geno.05.tmp

# 5% IID missingness threshold
plink \
  --bfile data/tmp/"$prefix"_qc1_geno.05.tmp \
  --mind 0.05 \
  --make-bed \
  --out data/tmp/"$prefix"_qc1_mind.05.tmp

# edit: delayed HWE until prep for PC-AiR
# # filter variants by HWE test
# # include non-controls because plink phenotypes aren't accurate
# plink \
#   --bfile data/tmp/"$prefix"_qc1_mind.05.tmp \
#   --hwe 1e-6 midp include-nonctrl \
#   --make-bed \
#   --out data/tmp/"$prefix"_qc1.tmp

# create pruned SNP subset
# exclude regions of known high heterozygosity
plink \
  --bfile data/tmp/"$prefix"_qc1_mind.05.tmp \
  --exclude raw_data/inversion.regions.txt \
  --indep-pairwise 15000 1500 0.2 \
  --out data/tmp/"$prefix"_qc1_pruned.tmp

# create list of participants with unusually high heterozygosity
plink \
  --bfile data/tmp/"$prefix"_qc1_mind.05.tmp \
  --extract data/tmp/"$prefix"_qc1_pruned.tmp.prune.in \
  --het \
  --out data/tmp/"$prefix"_qc1_het.tmp

# check heterozygosity and create list of participants who fail QC
# (absolute heterozygosity rate > 3 sd from mean)
Rscript --vanilla code/prep/act/03b_check_heterozygosity.R \
  data/tmp/"$prefix"_qc1_het.tmp.het \
  data/tmp/"$prefix"_qc1_fail_het.tmp

# create PLINK fileset that keeps only those who pass heterozygosity QC
plink \
  --bfile data/tmp/"$prefix"_qc1_mind.05.tmp \
  --remove ./data/tmp/"$prefix"_qc1_fail_het.tmp \
  --make-bed \
  --out data/tmp/"$prefix"_qc1.tmp

