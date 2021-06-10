#!/bin/bash

#==================================================
# list NACCIDs in ADGC imputed using TOPMed panel
# (run in gazelle)
#==================================================

for i in {1..12}
  do 
  file=$(ls /data_global/ADGC_GWAS/ADGC_NHW/ADC$i/RawGenotypes/ | egrep '.*.fam')
  cat /data_global/ADGC_GWAS/ADGC_NHW/ADC$i/RawGenotypes/$file | awk '{print$1,$2}' >> data/nacc_ids_adgc_topmed.txt
  done