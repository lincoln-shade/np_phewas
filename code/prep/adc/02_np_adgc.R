
library(data.table)
library(magrittr)

na_strings <- c("-9", "-4", "9999", "999", "888")
# Universal Data Set + NP
uds <- fread("/data_global/nacc/investigator_nacc56.csv", 
             na.strings = na_strings)
# NP data not longitudinal, so just keep obs from last visit
uds <- uds[NACCVNUM == NACCAVST]

# Minimal Data Set + NP
mds <- fread("/data_global/nacc/fardo09062019.csv", 
             na.strings = na_strings)

#------------------------------------------------------
# subset participants with NP data and NP variables
#------------------------------------------------------

# all autopsied have NPFORMVER
uds[!(is.na(NPFORMVER)), np := 1]
table(uds$np, useNA = "a")

mds[!(is.na(NPFORMVER)), np := 1]
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
np_mds_vars <- np_vars[!(np_vars %in% 
                           c("NPARTAG", "NPATGSEV", "NPATGAMY",
                             "NPATGAM1", "NPATGAM2", "NPATGAM3", 
                             "NPATGAM4", "NPATGAM5", "NPATGFRN", 
                             "NPATGFR1", "NPATGFR2", "NPATGFR3", 
                             "NPATGFR4"))]
mds_np <- mds[np == 1, ..np_mds_vars]
mds_np[, uds := 0]

np <- rbind(uds_np, mds_np, fill=TRUE)
# one duplicate ID, same values in MDS and UDS so I'll keep UDS row
sum(duplicated(np$NACCID))
np_dup_id <- np[NACCID %in% NACCID[duplicated(NACCID)], 
                -c("NACCMOD", "NACCYOD", "uds")]
np <- np[!((NACCID %in% np_dup_id$NACCID) & uds == 0)]

#-------------------------------------------------------------
# add variable indicating genotype data available in ADGC_HRC
#-------------------------------------------------------------

adgc <- fread("tmp/adc_no_naccid_dups.fam", header = FALSE)
np_adgc <-  merge(np, adgc[, .(V1, V2)], by.x = "NACCID", by.y = "V2")
np_adgc <- np_adgc[!(NACCID %in% NACCID[duplicated(NACCID)])]
setcolorder(np_adgc, c("V1", "NACCID"))
setnames(np_adgc, c("V1", "NACCID"), c("FID", "IID"))
saveRDS(np_adgc, file = "data/adc/np.Rds")
fwrite(np_adgc, file = "data/adc/np.csv", quote = FALSE)

