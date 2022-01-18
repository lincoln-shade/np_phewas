#!/bin/bash

#====================
# find related pairs
#====================
plink --bfile data/adc/adc_np_pruned --genome --min 0.18 --out data/adc/adc_np