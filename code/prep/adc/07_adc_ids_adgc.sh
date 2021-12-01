#!/bin/bash

#==================================================
# list NACCIDs in ADGC imputed using TOPMed panel
# (run in gazelle)
#==================================================

for i in {1..12}
  do 
  file=$(ls /data_global/ADGC_GWAS/ADGC_NHW/ADC$i/RawGenotypes/ | egrep '.*.fam')
  cat /data_global/ADGC_GWAS/ADGC_NHW/ADC$i/RawGenotypes/$file | awk -v adc=$i '{print$1,$2,adc}' >> data/adc_ids_adgc.txt
  done