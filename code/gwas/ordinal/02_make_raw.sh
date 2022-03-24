#!/bin/bash

<data/tmp/snp_list_index.tmp xargs -n 1 -P 22 -I % \
  plink \
    --bfile "$1" \
    --extract data/tmp/ordinal_snp_list_%.tmp \
    --allow-no-sex \
    --recode A \
    --out data/tmp/ordinal_snp_list_%
    