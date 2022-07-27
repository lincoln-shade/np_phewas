#!/bin/bash

pheno=$1
covar=$2
gene_annot=$3
out=$4


magma \
  --bfile data/mega/mega_np \
  --pheno file="$pheno" \
  --covar file="$pheno" \
  --gene-annot "$gene_annot" \
  --out "$out"