#!/bin/bash

##--------------------------------------------
## generate list of independent top variants
##--------------------------------------------



plink \
  --bfile $1 \
  --clump $2 \
  --clump-p1 0.00001 \
  --clump-r2 0.05 \
  --out $3