#!/bin/bash

#==================================================
# list IIDs in ADGC imputed using TOPMed panel
# (run in gazelle)
#==================================================
rm data/act/act_ids_adgc.txt
for i in {1..3}
  do 
  file=$(ls /data_global/ADGC_GWAS/ADGC_NHW/ACT$i/CleanedGenotypes/ | egrep '.*.fam')
  cat /data_global/ADGC_GWAS/ADGC_NHW/ACT$i/CleanedGenotypes/$file | \
  awk -v cohort=$i '{print$1,$2,cohort}' >> data/act/act_ids_adgc.txt
  done