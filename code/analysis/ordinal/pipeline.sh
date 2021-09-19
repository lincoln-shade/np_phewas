#!/bin/bash

cat data/adc_np_ord.txt | awk 'NR==1{for (i=3; i<=NF; i++) print $i}' | xargs -n 1 -P 4 bash code/analysis/ordinal/03_regression.sh