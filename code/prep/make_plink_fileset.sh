#!/bin/bash

for i in {1..20}
  do
  bash code/prep/vcf_to_plink.sh $i
  done

# rm data/tmp/ADC[0-9]*