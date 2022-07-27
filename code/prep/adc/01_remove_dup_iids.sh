#!/bin/bash

#==========================================
# basic QC
#==========================================

pfile=tmp/adc_rsid
prefix=adc

#--------------------------------
# identify duplicates
#--------------------------------

plink2 --pfile "$pfile" --maf 0.05 --make-bed --out tmp/"$prefix"_maf05

king -b tmp/"$prefix"_maf05.bed --duplicate --prefix tmp/"$prefix"_dup

cat tmp/"$prefix"_dup.con | \
    awk '{if ($2 != "ID1") {print $1, $2} }' \
    > tmp/"$prefix"_dup_id.tmp

cat tmp/"$prefix"_dup.con | \
    awk '{if ($4 != "ID2") {print $3, $4} }' \
    >> tmp/"$prefix"_dup_id.tmp

plink \
    --bfile tmp/"$prefix"_maf05 \
    --missing \
    --keep tmp/"$prefix"_dup_id.tmp \
    --out tmp/"$prefix"_dup

Rscript --vanilla code/prep/adc/01a_dup_list.R

plink2 \
    --pfile tmp/"$prefix"_rsid \
    --remove tmp/"$prefix"_dup_remove.tmp \
    --make-pgen \
    --out tmp/"$prefix"_no_dup

# now need to remove duplicant NACCIDs from .psam
Rscript --vanilla code/prep/adc/01b_dup_naccid_list.R

naccid_dups=tmp/adc_dup_naccid.tmp
pfile_no_naccid_dups=tmp/"$prefix"_no_naccid_dups
if test -f  "$naccid_dups"; 
then
    plink2 \
        --pfile tmp/"$prefix"_no_dup \
        --remove "$naccid_dups" \
        --make-bed \
        --out "$pfile_no_naccid_dups"
else
    mv tmp/"$prefix"_no_dup.pgen "$pfile_no_naccid_dups".pgen
    mv tmp/"$prefix"_no_dup.pvar "$pfile_no_naccid_dups".pvar
    mv tmp/"$prefix"_no_dup.psam "$pfile_no_naccid_dups".psam
    
fi

# rename .fam IIDs to NACCIDs
Rscript --vanilla code/prep/adc/01c_rename_psam.R

mv "$pfile_no_naccid_dups".pgen data/adc/adc.pgen
mv "$pfile_no_naccid_dups".pvar data/adc/adc.pvar
mv "$pfile_no_naccid_dups".psam data/adc/adc.psam

