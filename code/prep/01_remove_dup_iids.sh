#!/bin/bash

#==========================================
# basic QC
#==========================================

bfile=data/plink/adc

#--------------------------------
# identify duplicates
#--------------------------------

plink --bfile "$bfile" --maf 0.05 --make-bed --out data/tmp/adc_maf05

king -b data/tmp/adc_maf05.bed --duplicate --prefix data/tmp/adc_dup

cat data/tmp/adc_dup.con | awk '{if ($2 != "ID1") {print $1, $2} }' >> data/tmp/adc_dup_id.tmp
cat data/tmp/adc_dup.con | awk '{if ($4 != "ID2") {print $3, $4} }' >> data/tmp/adc_dup_id.tmp

plink --bfile data/tmp/adc_maf05 --missing --keep data/tmp/adc_dup_id.tmp --out data/tmp/adc_dup

Rscript --vanilla code/prep/01a_make_dup_id_list.R

plink --bfile data/plink/adc --remove data/tmp/adc_dup_remove.tmp --make-bed --out data/tmp/adc_no_dup
# now need to remove "_dup" from iids in R
Rscript --vanilla code/prep/01b_remove_id_dup_suffix_fam.R
