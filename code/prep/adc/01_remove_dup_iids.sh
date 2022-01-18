#!/bin/bash

#==========================================
# basic QC
#==========================================

bfile=data/adc/adc
prefix=adc

#--------------------------------
# identify duplicates
#--------------------------------

plink --bfile "$bfile" --maf 0.05 --make-bed --out tmp/"$prefix"_maf05

king -b tmp/"$prefix"_maf05.bed --duplicate --prefix tmp/"$prefix"_dup

cat tmp/"$prefix"_dup.con | awk '{if ($2 != "ID1") {print $1, $2} }' >> tmp/"$prefix"_dup_id.tmp
cat tmp/"$prefix"_dup.con | awk '{if ($4 != "ID2") {print $3, $4} }' >> tmp/"$prefix"_dup_id.tmp

plink --bfile tmp/"$prefix"_maf05 --missing --keep tmp/"$prefix"_dup_id.tmp --out tmp/"$prefix"_dup

Rscript --vanilla code/prep/adc/01a_dup_list.R

plink --bfile data/plink/"$prefix" --remove tmp/"$prefix"_dup_remove.tmp --make-bed --out tmp/"$prefix"_no_dup
# now need to remove NACCIDs from .fam
Rscript --vanilla code/prep/adc/01b_dup_naccid_list.R

naccid_dups=tmp/adc_dup_naccid.tmp
bfile_no_naccid_dups=tmp/adc_no_naccid_dups
if test -f  "$naccid_dups"; 
then
    plink --bfile tmp/"$prefix"_no_dup --remove "$naccid_dups" --make-bed --out "$bfile_no_naccid_dups"
else
    mv tmp/"$prefix"_no_dup.bed "$bfile_no_naccid_dups".bed
    mv tmp/"$prefix"_no_dup.bim "$bfile_no_naccid_dups".bim
    mv tmp/"$prefix"_no_dup.fam "$bfile_no_naccid_dups".fam
    
fi

# rename .fam IIDs to NACCIDs
Rscript --vanilla code/prep/adc/01c_rename_fam.R

