#!/bin/bash

# seq 1 22 | xargs -n 1 -P 22 -I % plink --vcf /data_global/ADC_imputed/chr%_all_adcs.dose.vcf.gz --make-bed --const-fid --out data/tmp/adc_chr%

# seq 1 22 | xargs -n 1 -P 22 -I % plink --bfile data/tmp/adc_chr% --maf 0.001 --geno 0.2 --make-bed --out data/tmp/adc_chr%_qc

# list_file=data/tmp/plink_prefix_list.tmp
# echo > $list_file
# for i in $(seq 1 22)
#   do
#   echo data/tmp/adc_chr"$i"_qc >> $list_file
#   done
# 
# plink --merge-list $list_file --make-bed --out data/plink/adc

plink --bfile data/plink/adc --update-name raw_data/chrall.rs.linker --make-bed --out data/plink/adc_rsid