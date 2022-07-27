#===============================================
# Add systolic blood pressure to covariate file
#===============================================

library(data.table)
library(readxl)

covar <- fread("data/mega/mega_np_age80.covar")

na_strings <- c("-9", "-4", "9999", "999", "888")
uds <- fread("/data_global/nacc/investigator_nacc57.csv", 
             na.strings = na_strings)
# NP data not longitudinal, so just keep obs from last visit
uds <- uds[NACCVNUM == NACCAVST]
uds <- uds[!is.na(BPSYS), .(NACCID, BPSYS)]
setnames(uds, colnames(uds), c("IID", "sbp"))

rosmap <- setDT(read_xlsx(paste0("/data_global/ROSMAP/greg_20200109/",
                                 "dataset_843_long_01-09-2020.xlsx")))
most_recent_visit_ind <- 
  rosmap[, .I[age_at_visit == max(age_at_visit)], projid]$V1

rosmap <- rosmap[most_recent_visit_ind, .(projid, sbp_avg)]
rosmap <- rosmap[!is.na(sbp_avg)]
setnames(rosmap, colnames(rosmap), c("IID", "sbp"))

sbp <- rbind(uds, rosmap)

covar_sbp <- merge(covar, sbp, by = "IID")
# hypotension SBP 90 mmHg or below, remove so that association
# more likely to be monotonic
covar_sbp <- covar_sbp[sbp >= 90]

setcolorder(covar_sbp, c("FID", "IID"))
fwrite(covar_sbp, file = "data/mega/mega_np_age80_sbp.covar", sep = " ", 
       quote = FALSE)

