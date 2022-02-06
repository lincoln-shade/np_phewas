#====================================================
# create list of ADNI participants with autopsy data
#====================================================

library(data.table)
library(magrittr)
library(stringi)

np <- fread("raw_data/NEUROPATH_07_06_21.csv")
fam <- fread("/data_global/ADNI/ADNI_TOPMed_fromHohman_20211214/ADNI_NHW_imputed_final.fam")
ids <- merge(fam, np, by.x = 'V2', by.y = 'RID')
ids[, RID := as.integer(RID)]
