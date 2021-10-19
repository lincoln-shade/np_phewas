!#/bin/bash

pheno1=$1
pheno2=$2
out=output/gen_cor/"$pheno1"_"$pheno2"
summary1=data/sumstats/"$pheno1"_sumstats.txt
summary2=data/sumstats/"$pheno2"_sumstats.txt
ldak5.1.linux --sum-cors $out --summary $summary1 --summary2 $summary2 --tagfile raw_data/ldak.thin.hapmap.gbr.tagging --allow-ambiguous YES --cutoff 0.01 --check-sums NO