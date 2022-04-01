#!/bin/bash

cat data/mega/mega_np.pheno | awk 'NR==1{for (i=3; i<=NF; i++) print $i}' | xargs -n 1 bash code/coloc/00_coloc_pipeline.sh