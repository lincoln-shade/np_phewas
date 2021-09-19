#!/bin/bash

#------------------------------
# ROSMAP ordinal regression
#------------------------------

# Rscript code/single_variant_analysis/ordinal/01_snp_list.R data/rosmap/rosmap_unrelated.bim

# bash code/single_variant_analysis/ordinal/02_regression.sh data/rosmap/rosmap_unrelated.RData "arteriol_scler" data/rosmap/rosmap_unrelated

Rscript code/single_variant_analysis/ordinal/04_results.R "output/rosmap/rosmap.assoc.logistic" "output/rosmap/"