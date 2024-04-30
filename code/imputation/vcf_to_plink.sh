#!/bin/bash

# # convert imputed VCF files to PLINK 2 binary filesets
# seq 1 22 | \
#   xargs -n 1 -P 22 -I % \
#   ./code/imputation/vcf_to_plink2.py \
#     --vcf /data_global/ADC_imputed/chr%_all_adcs.dose.vcf.gz \
#     --out tmp/adc_chr%

# merge PLINK 2 chromosomal binary filesets
list_file=tmp/adc_prefix_list.tmp
# for i in $(seq 1 22)
#   do
#   echo tmp/adc_chr"$i" >> $list_file
#   done

plink2 --pmerge-list $list_file --out tmp/adc

# rename variants with rsIDs (doesn't work for indels afaict)
linker=tmp/rsid_linker.tmp
Rscript --vanilla ./code/imputation/import_rsids_from_bed.R \
  --pvar tmp/adc.pvar \
  --out "$linker"

plink2 --pfile tmp/adc --update-name "$linker" --make-pgen --out tmp/adc_rsid