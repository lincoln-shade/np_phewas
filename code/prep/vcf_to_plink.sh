#!/bin/bash

#=========================================================================
# create plink filesets for each autosomal chromosome in folder,
# then merge into one plink fileset
#=========================================================================

# note, not all combination of cohort/race present
study=ADC
chr=$1 # 1, 2, or 3
race=NHW # NHW, Asian, Hispanic, or AA

# make folder if needed
if [[ ! -d "data" ]]
  then
  mkdir data
fi

if [[ ! -d "data/tmp" ]]
  then
  mkdir data/tmp
fi

vcf_prefix_1=/data_global/ADGC_GWAS/ADGC_"$race"/"$study"

if [[ $race = "NHW" ]]
then 
  vcf_prefix_2=/TOPMEDr2/vcffile/chr"$chr".dose.vcf.gz
else 
  vcf_prefix=_"$race"/TOPMEDr2/vcffile/chr"$chr".dose.vcf.gz
fi

# create file to list fileset prefixes to merge
list_file=data/tmp/"$study"_"$race"_chr"$chr"_plink_fileset_prefixes.tmp.txt
echo > $list_file
seq 1 12 | xargs -I % echo ./data/tmp/"$study"%_"$race"_chr"$chr" >> $list_file

# create plink fileset for each chromosome.
seq 1 12 | xargs -P 12 -I %  python3 code/prep/convert_vcf.py --vcf "$vcf_prefix_1"%"$vcf_prefix_2" --geno 0.5 --out ./data/tmp/"$study"%_"$race"_chr"$chr"


plink \
  --merge-list $list_file \
  --make-bed \
  --maf 0.001 \
  --geno 0.2 \
  --out data/tmp/"$study"_"$race"_chr"$chr"
