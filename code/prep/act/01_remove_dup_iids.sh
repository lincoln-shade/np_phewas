#!/bin/bash

#==========================================
# basic QC
#==========================================

bfile=data/act/act
prefix=act

#--------------------------------
# identify duplicates
#--------------------------------

plink --bfile "$bfile" --maf 0.05 --make-bed --out tmp/"$prefix"_maf05

king -b tmp/"$prefix"_maf05.bed --duplicate --prefix tmp/"$prefix"_dup

cat tmp/"$prefix"_dup.con | awk '{if ($2 != "ID1") {print $1, $2} }' >> tmp/"$prefix"_dup_id.tmp
cat tmp/"$prefix"_dup.con | awk '{if ($4 != "ID2") {print $3, $4} }' >> tmp/"$prefix"_dup_id.tmp

plink --bfile tmp/"$prefix"_maf05 --missing --keep tmp/"$prefix"_dup_id.tmp --out tmp/"$prefix"_dup

Rscript --vanilla code/prep/act/01a_make_dup_id_list.R $prefix

plink --bfile data/act/"$prefix" --remove tmp/"$prefix"_dup_remove.tmp --make-bed --out tmp/"$prefix"_no_dup
# now need to remove from iids in R
Rscript --vanilla code/prep/act/01b_remove_id_dup_suffix_fam.R

# rename SNPs with rsIDs
Rscript --vanilla code/prep/act/01c_update_rsid_linker.R

plink --bfile tmp/"$prefix"_no_dup --update-name tmp/rsid_linker.tmp --make-bed --out tmp/"$prefix"_no_dup_rsid