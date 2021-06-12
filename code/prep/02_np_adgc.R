
library(pacman)
p_load("data.table", "magrittr", "ggplot2", "pacman")
na_strings <- c("-9", "-4", "9999", "999", "888")
# Universal Data Set + NP
uds <- fread("/data_global/nacc/investigator_nacc53.csv", na.strings = na_strings)
# NP data not longitudinal, so just keep obs from last visit
uds <- uds[NACCVNUM == NACCAVST]


# Minimal Data Set + NP
mds <- fread("/data_global/nacc/fardo09062019.csv", na.strings = na_strings)
#------------------------------------------------------
# subset participants with NP data and NP variables
#------------------------------------------------------

# From the NACC NP data dictionary, regarding NACCDAGE:
# > "Note that all subjects with an NP Form will have a known age at death."
# So, can use NACCDAGE to separate those with NP data from those without

uds[!(is.na(NACCDAGE)), np := 1]
table(uds$np, useNA = "a")

mds[!(is.na(NACCDAGE)), np := 1]
table(mds$np, useNA = "a")

np_vars <- colnames(uds)[grep("NP", colnames(uds))]
np_vars_nacc <- fread("data/nacc_np_vars_names_nacc.txt", header = FALSE)$V1
np_vars_remove <- fread("data/nacc_uds_var_names_np.txt", header = FALSE)$V1
np_vars <- np_vars[!(np_vars) %in% np_vars_remove]
np_vars <- c(np_vars_nacc, np_vars)
uds_np <- uds[np == 1, ..np_vars]
uds_np[, uds := 1]
# month of death and year of death not in MDS, but still want to have for UDS
mds[, `:=`(NACCMOD = NA, NACCYOD = NA)]
mds_np <- mds[np == 1, ..np_vars]
mds_np[, uds := 0]

np <- rbind(uds_np, mds_np)
# one duplicate ID, same values in MDS and UDS so I'll keep UDS row
sum(duplicated(np$NACCID))
np_dup_id <- np[NACCID %in% NACCID[duplicated(NACCID)], -c("NACCMOD", "NACCYOD", "uds")]
np <- np[!((NACCID %in% np_dup_id$NACCID) & uds == 0)]

#-------------------------------------------------------------
# add variable indicating genotype data available in ADGC_HRC
#-------------------------------------------------------------

adgc_hrc <- fread("/data_global/ADGC_HRC/converted/full/adgc_hrc_merged_qced.covar")
adgc_raw <- fread("data/nacc_ids_adgc_topmed.txt", header = FALSE)
adgc_topmed <- fread("data/ADC_NHW_chr1.fam", header = FALSE)
np_hrc <- np[NACCID %in% adgc_hrc$IID]
np_topmed <-  merge(np, adgc_topmed[, .(V1, V2)], by.x = "NACCID", by.y = "V2")
np_topmed <- np_topmed[!(NACCID %in% NACCID[duplicated(NACCID)])]
np_hrc_only <- np_hrc[!(NACCID %in% np_topmed$NACCID)]
np_topmed_only <- np_topmed[!(NACCID %in% np_hrc$NACCID)]
setcolorder(np_topmed, c("V1", "NACCID"))
setnames(np_topmed, c("V1", "NACCID"), c("FID", "IID"))
saveRDS(np_topmed, file = "data/nacc_np.Rds")
rm(list = ls())
p_unload(all)
