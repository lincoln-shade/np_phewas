#!/bin/bash

#==========================================
# check relatedness in NACC/ADGC dataset
# and perform other QC measures
#==========================================


# create PLINK fileset of participants with B-ASC phenotyping
plink \
  --bfile /data_global/ADGC_HRC/converted/full/adgc_hrc_merged_qced \
  --keep data/tmp/nacc_ids.tmp.txt \
  --make-bed \
  --out data/tmp/nacc.tmp

#-------------------------------------------------------------
# Step 1: basic variant and individual quality filters
# remove participants with missing genotyping rate of > 0.05
# remove variants with MAF < 0.05
# remove variants that fail HWE test
# remove variants that have missing genotype rate > 0.05
#--------------------------------------------------------------

# 20% variant missingness threshold
plink \
  --bfile data/tmp/nacc.tmp \
  --geno 0.2 \
  --make-bed \
  --out data/tmp/nacc_qc1_geno.2.tmp

# 20% IID missingness threshold
plink \
  --bfile data/tmp/nacc_qc1_geno.2.tmp \
  --mind 0.2 \
  --make-bed \
  --out data/tmp/nacc_qc1_mind.2.tmp

# 5% variant missingness threshold
plink \
  --bfile data/tmp/nacc_qc1_mind.2.tmp \
  --geno 0.05 \
  --make-bed \
  --out data/tmp/nacc_qc1_geno.05.tmp

# 5% IID missingness threshold
plink \
  --bfile data/tmp/nacc_qc1_geno.05.tmp \
  --mind 0.05 \
  --make-bed \
  --out data/tmp/nacc_qc1_mind.05.tmp

# filter variants by HWE test
# include non-controls because plink phenotypes aren't accurate
plink \
  --bfile data/tmp/nacc_qc1_mind.05.tmp \
  --hwe 1e-6 midp include-nonctrl \
  --make-bed \
  --out data/tmp/nacc_qc1.tmp

# create pruned SNP subset
# exclude regions of known high heterozygosity
plink \
  --bfile data/tmp/nacc_qc1.tmp \
  --exclude data/inversion.regions.txt \
  --indep-pairwise 15000 1500 0.2 \
  --out data/tmp/nacc_qc1_pruned.tmp

# create list of participants with unusually high heterozygosity
plink \
  --bfile data/tmp/nacc_qc1.tmp \
  --extract data/tmp/nacc_qc1_pruned.tmp.prune.in \
  --het \
  --out data/tmp/nacc_qc1_het.tmp

# check heterozygosity and create list of participants who fail QC
# (absolute heterozygosity rate > 3 sd from mean)
Rscript code/prep/nacc/02a_check_heterozygosity.R

# create PLINK fileset that keeps only those who pass heterozygosity QC
plink \
  --bfile data/tmp/nacc_qc1.tmp \
  --remove ./data/tmp/nacc_qc1_fail_het.tmp.txt \
  --make-bed \
  --out data/tmp/nacc_qc2.tmp

