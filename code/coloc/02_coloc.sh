#!/bin/bash

#====================================
# Run colocalization analysis
#====================================

# command-line arguments
# 1: PLINK binary fileset prefix
# 2: PLINK .pheno file name
# 3: PLINK .assoc.logistic file name
# 4: relative path for colocalization results output directory 
# (append with "/" e.g. ./ or data/)

# calculate minor allele frequencies for each chromosome 
# that has at least 1 qtl
pfile=$1
pheno=$2
regression=$3
out_folder=$4
gwas_pheno=$5

chrs=tmp/chrs_"$gwas_pheno".tmp
maf_out=tmp/chr%_maf_"$gwas_pheno"
cat $chrs | xargs -n 1 -P 10 -I % \
  plink2 \
    --pfile $pfile \
    --chr % \
    --freq \
    --out "$maf_out"

# prepare gwas summary stats data


cat $chrs | xargs -n 1 -P 10 -I % \
  Rscript code/coloc/02a_gwas_sumstats.R \
    -p $pheno \
    -f "$maf_out".afreq \
    -r $regression \
    --chr % \
    --phenotype $gwas_pheno

# merge gwas and qtl data
cat tmp/top_qtls_"$gwas_pheno".tmp | \
  awk 'NR > 1 {print$1,$2,$3,$4}' | \
  xargs -P 22 -L 1 bash -c \
  'code/coloc/02b_gwas_qtl.R -g $0 -q $1 -t $2 --chr $3 --qtl $4' $gwas_pheno


# run coloc
p12=0.00001
cat tmp/top_qtls_"$gwas_pheno".tmp | \
  awk 'NR>1 {print$1,$2,$3}' | \
  xargs -L 1 bash -c \
  './code/coloc/02c_coloc.R -f tmp/chr"$5"_"$3"_"$4"_gwas_qtl_"$2".tmp --p12 $0 -o "$1""$2"_coloc_results.txt -q $3 -t $4 --chr $5' \
  $p12 $out_folder $gwas_pheno
