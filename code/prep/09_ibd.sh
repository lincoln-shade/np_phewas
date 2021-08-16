#!/bin/bash

#====================
# find related pairs
#====================
plink --bfile data/plink/adc_np_pruned --genome --min 0.18 --out data/plink/adc_np