#!/bin/bash
out_prefix=output/act/act_mega
data_prefix=data/act/act_mega
# bash code/coloc/01_qtls.sh "$out_prefix".assoc.logistic
bash code/coloc/02_coloc.sh $data_prefix "$data_prefix".pheno "$out_prefix".assoc.logistic output/act/