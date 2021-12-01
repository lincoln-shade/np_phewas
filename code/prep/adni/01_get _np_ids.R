#====================================================
# create list of ADNI participants with autopsy data
#====================================================

library(data.table)
library(magrittr)
library(stringi)

np <- fread("raw_data/NEUROPATH_07_06_21.csv")

fam <- fread("raw_data/WGS_Omni25_BIN_wo_ConsentsIssues.fam")
adni1_fam <- fread("raw_data/ADNI_cluster_01_forward_757LONI.fam")
adni2_fam <- fread("raw_data/ADNI_GO_2_Forward_Bin.fam")
adni3_fam <- fread("raw_data/PLINK_Final/ADNI3_PLINK_Final.fam")
roster <- fread("raw_data/ROSTER.csv")

keep_ids <- roster[RID %in% np$RID, .(RID, PTID)] %>% .[!duplicated(.)]

merge_ids_fam <- function(dt) {
  merge(dt[, .(V1, V2)], keep_ids, by.x = c("V2"), by.y = c("PTID"))
}

ids <- rbindlist(list(merge_ids_fam(fam), 
                      merge_ids_fam(adni1_fam),
                      merge_ids_fam(adni2_fam),
                      merge_ids_fam(adni3_fam)))

ids[, RID := as.integer(RID)]
