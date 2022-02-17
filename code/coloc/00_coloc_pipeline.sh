#!/bin/bash
phenotype=$1
out_prefix=output/gwas/mega/"$phenotype"
data_prefix=data/mega/mega_np
bash code/coloc/01_qtls.sh $phenotype $out_prefix

top_qtls_file=data/tmp/top_qtls_"$phenotype".tmp

if test -f $top_qtls_file; then
    bash code/coloc/02_coloc.sh $data_prefix "$data_prefix".pheno "$out_prefix".assoc.logistic output/coloc/mega/"$phenotype" $phenotype
fi