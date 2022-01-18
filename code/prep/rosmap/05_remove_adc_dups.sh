#!/bin/bash

# remove ROSMAP participants that are duplicates of ADC

plink \
  --bfile tmp/rosmap_2 \
  --keep data/rosmap/ids_qced.txt \
  --make-bed \
  --out tmp/rosmap_3

# merge adc and 1000g filesets
python code/prep/merge_plink_filesets.py -b tmp/rosmap_3 -b1 data/adc/adc_np -o tmp/rosmap_adc_merged

king -b tmp/rosmap_adc_merged.bed --duplicate --prefix tmp/rosmap_adc_dup

# no duplicates found

# check relatedness

king -b tmp/rosmap_adc_merged.bed --ibdseg --degree 2 --prefix data/rosmap/rosmap_adc