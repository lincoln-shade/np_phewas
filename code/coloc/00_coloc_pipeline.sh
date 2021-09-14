#!/bin/bash
phenotype=$1
out_prefix=output/gwas_results/adc_np_"$phenotype"
data_prefix=data/plink/adc_np
bash code/coloc/01_qtls.sh $phenotype

top_qtls_file=data/tmp/top_qtls_"$phenotype".tmp

if test -f $top_qtls_file; then
    bash code/coloc/02_coloc.sh $data_prefix "$data_prefix".pheno "$out_prefix".assoc.logistic output/coloc/"$phenotype" $phenotype
fi