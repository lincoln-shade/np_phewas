#!/bin/bash

# for i in {1..22}
#   do
#   bash code/prep/vcf_to_plink.sh $i
#   rm data/tmp/ADC[0-9]*chr$i.*
#   done

list_file=data/tmp/ADC_NHW_plink_fileset_prefixes.tmp
echo > $list_file
for i in {1..22}
  do
  echo data/tmp/ADC_NHW_chr$i >> $list_file
  done

if [[ ! -d "data/plink" ]]
  then
  mkdir data/plink
fi

plink \
  --merge-list $list_file \
  --make-bed \
  --maf 0.001 \
  --geno 0.2 \
  --out data/plink/nacc