#!/bin/bash

#------------------------------
# NACC/ADGC ordinal regression
#------------------------------

Rscript code/single_variant_analysis/ordinal/01_snp_list.R data/nacc_adgc/nacc_adgc_unrelated.bim 

bash code/single_variant_analysis/ordinal/02_regression.sh data/nacc_adgc/nacc_adgc_unrelated.RData "NACCARTE"