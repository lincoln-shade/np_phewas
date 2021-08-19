#!/bin/bash

#====================================
# Run colocalization analysis
#====================================

# command-line arguments
# 1: PLINK binary fileset prefix
# 2: PLINK .pheno file name
# 3: PLINK .assoc.logistic file name
# 4: relative path for colocalization results output directory (append with "/" e.g. ./ or data/)

# calculate minor allele frequencies for each chromosome that has at least 1 qtl
bfile=$1
chrs=data/tmp/chrs.tmp
cat $chrs | xargs -n 1 -P 10 -I % plink --bfile $bfile --chr % --freq --out "data/tmp/chr%_maf"

# prepare gwas summary stats data
pheno=$2
regression=$3
cat $chrs | xargs -n 1 -P 10 -I % Rscript code/coloc/02a_gwas_sumstats.R $pheno "data/tmp/chr%_maf.frq" $regression %

# merge gwas and qtl data
cat data/tmp/top_qtls.tmp | awk 'NR>1 {print$1,$2,$3,$4}' | xargs -P 10 -L 1 Rscript code/coloc/02b_gwas_qtl.R

# run coloc
out_folder=$4
p12=0.00001
cat data/tmp/top_qtls.tmp | awk 'NR>1 {print$1,$2,$3}' | xargs -L 1 Rscript code/coloc/02c_coloc.R $p12 $out_folder
