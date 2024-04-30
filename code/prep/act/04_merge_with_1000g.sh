#!/bin/bash

#================================
# Step 3: Remove ethnic outliers
#================================
prefix=act
#---------------------------------------
# Sync variant IDs between filesets
#---------------------------------------

# need to remove :A1:A2 from "$prefix" variant IDs
cat tmp/"$prefix"_qc1.tmp.bim | awk 'gsub(/:[ACTG]*:[ACTG]*/, "", $2)' > tmp/"$prefix"_qc1_varID_noA1A2.tmp.bim

# write duplicate variant ids to file
cat tmp/"$prefix"_qc1_varID_noA1A2.tmp.bim | cut -d ' ' -f 2 | sort | uniq -d > tmp/"$prefix"_qc1_dups.tmp.txt

# create plink fileset with only biallelic variants
plink \
  --bfile tmp/"$prefix"_qc1.tmp \
  --bim tmp/"$prefix"_qc1_varID_noA1A2.tmp.bim \
  --exclude tmp/"$prefix"_qc1_dups.tmp.txt \
  --make-bed \
  --out tmp/"$prefix"_qc2a.tmp

# create list of SNPs from QCed "$prefix"
awk '{print$2}' tmp/"$prefix"_qc2a.tmp.bim > tmp/"$prefix"_qc2a_vars.tmp.txt

# create 1000g fileset with "$prefix" snps
plink \
  --bfile data/1000g/1000g \
  --extract tmp/"$prefix"_qc2a_vars.tmp.txt \
  --make-bed \
  --out tmp/1000g1.tmp

# create list of filtered 1000g snps
awk '{print$2}' tmp/1000g1.tmp.bim > tmp/1000g1_vars.tmp.txt

# create "$prefix" fileset with filtered 1000g snps
plink \
  --bfile tmp/"$prefix"_qc2a.tmp \
  --biallelic-only strict \
  --extract tmp/1000g1_vars.tmp.txt \
  --make-bed \
  --out tmp/"$prefix"_qc1c.tmp

# now both filesets have the exact same variant ids

#------------------------------------------
# Ensure both filesets have same build
#------------------------------------------
# sync BPs between filesets
plink \
  --bfile tmp/1000g1.tmp \
  --update-map tmp/"$prefix"_qc1c.tmp.bim 4 2 \
  --make-bed \
  --out tmp/1000g2.tmp

#---------------------------------------------
# ensure all variants have same A1/A2 alleles
#---------------------------------------------

# 1 set reference genome to "$prefix" alleles
plink \
  --bfile tmp/1000g2.tmp \
  --a1-allele tmp/"$prefix"_qc1c.tmp.bim 5 2 \
  --make-bed \
  --out tmp/1000g3.tmp

# 2 try to resolve strand issues
awk '{print$2,$5,$6}' tmp/"$prefix"_qc1c.tmp.bim > tmp/"$prefix"_qc1c_var.tmp.txt
awk '{print$2,$5,$6}' tmp/1000g3.tmp.bim > tmp/1000g3_var.tmp.txt
sort tmp/"$prefix"_qc1c_var.tmp.txt tmp/1000g3_var.tmp.txt |uniq -u > tmp/"$prefix"_qc1c_1000g3_var_diff.tmp.txt

awk '{print$1}' tmp/"$prefix"_qc1c_1000g3_var_diff.tmp.txt | sort -u > tmp/"$prefix"_1000g_flip_list.tmp.txt

plink \
  --bfile tmp/1000g3.tmp \
  --flip tmp/"$prefix"_1000g_flip_list.tmp.txt \
  --a1-allele tmp/"$prefix"_qc1c.tmp.bim 5 2 \
  --make-bed \
  --out tmp/1000g4.tmp

awk '{print$2,$5,$6}' tmp/1000g4.tmp.bim > tmp/1000g4_var.tmp.txt
sort tmp/"$prefix"_qc1c_var.tmp.txt tmp/1000g4_var.tmp.txt |uniq -u  > tmp/"$prefix"_1000g_uncorrected_var.tmp.txt

# none of the variants were fixed by flipping.

# 3 remove uncorrected variants
awk '{print$1}' tmp/"$prefix"_1000g_uncorrected_var.tmp.txt | sort -u > tmp/"$prefix"_1000g_exclude_var.tmp.txt

plink \
  --bfile tmp/"$prefix"_qc1c.tmp \
  --exclude tmp/"$prefix"_1000g_exclude_var.tmp.txt \
  --make-bed \
  --out tmp/"$prefix"_qc1d.tmp

plink \
  --bfile tmp/1000g4.tmp \
  --exclude tmp/"$prefix"_1000g_exclude_var.tmp.txt \
  --make-bed \
  --out tmp/1000g5.tmp

#---------------------------
# merge "$prefix" and 1000g
#---------------------------

# merge
plink \
  --bfile tmp/"$prefix"_qc1d.tmp \
  --bmerge tmp/1000g5.tmp \
  --allow-no-sex \
  --make-bed \
  --out tmp/"$prefix"_1000g_merged

# prune
plink \
  --bfile tmp/"$prefix"_1000g_merged \
  --indep-pairwise 15000 1500 0.2 \
  --out tmp/"$prefix"_1000g_merged_prune

# pca
plink \
  --bfile tmp/"$prefix"_1000g_merged \
  --no-pheno \
  --extract tmp/"$prefix"_1000g_merged_prune.prune.in \
  --pca 2 \
  --out tmp/"$prefix"_1000g_merged_pca
  