#!/bin/bash

#==============================================================
# basic QC
#==============================================================

#--------------------------------------------------------------
# make plink fileset of only those with autopsy data
#--------------------------------------------------------------

cat data/adc/np.csv | awk -F "," 'NR>1{print $1,$2}' > tmp/np_ids.tmp

# apply initial MAF filter of 1%
plink2 --pfile data/adc/adc --keep tmp/np_ids.tmp --maf 0.01 --make-pgen --out tmp/adc_np

#-------------------------------------------------------------
# Step 1: basic variant and individual quality filters
# remove participants with missing genotyping rate of > 0.05
# remove variants with MAF < 0.05
# remove variants that fail HWE test
# remove variants that have missing genotype rate > 0.05
#--------------------------------------------------------------

# 20% variant missingness threshold
plink2 \
  --pfile tmp/adc_np \
  --geno 0.2 \
  --make-pgen \
  --out tmp/adc_qc1_geno.2.tmp

# 20% IID missingness threshold
plink2 \
  --pfile tmp/adc_qc1_geno.2.tmp \
  --mind 0.2 \
  --make-pgen \
  --out tmp/adc_qc1_mind.2.tmp

# 5% variant missingness threshold
plink2 \
  --pfile tmp/adc_qc1_mind.2.tmp \
  --geno 0.05 \
  --make-pgen \
  --out tmp/adc_qc1_geno.05.tmp

# 5% IID missingness threshold
plink2 \
  --pfile tmp/adc_qc1_geno.05.tmp \
  --mind 0.05 \
  --make-pgen \
  --out tmp/adc_qc1_mind.05.tmp

# create pruned SNP subset
# exclude regions of known high heterozygosity
plink2 \
  --pfile tmp/adc_qc1_mind.05.tmp \
  --exclude data/inversion.regions.txt \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/adc_qc1_pruned.tmp

# create list of participants with unusually high heterozygosity
plink2 \
  --pfile tmp/adc_qc1_mind.05.tmp \
  --extract tmp/adc_qc1_pruned.tmp.prune.in \
  --het \
  --out tmp/adc_qc1_het.tmp

# check heterozygosity and create list of participants who fail QC
# (absolute heterozygosity rate > 3 sd from mean)
Rscript --vanilla code/prep/adc/03b_check_heterozygosity.R

# create PLINK fileset that keeps only those who pass heterozygosity QC
plink2 \
  --pfile tmp/adc_qc1_mind.05.tmp \
  --remove ./tmp/adc_qc1_fail_het.tmp \
  --make-pgen \
  --out tmp/adc_qc1.tmp

