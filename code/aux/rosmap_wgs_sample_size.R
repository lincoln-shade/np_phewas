library(data.table)
library(stringi)
library(readxl)

linker <- fread("~/ROSMAP_biospecimen_metadata.csv")
wgs <- fread("~/ROSMAP_assay_wholeGenomeSeq_metadata.csv")
clinical <- fread("~/ROSMAP_clinical.csv")
np <- setDT(read_xlsx(
  "/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))

linker_wgs <- merge(linker[assay == "wholeGenomeSeq", 
                           .(specimenID, individualID)], 
                    wgs[, .(specimenID)], 
                    by = "specimenID")
linker_wgs <- merge(linker_wgs[, .(individualID)], 
                    clinical[, .(individualID, projid)], "individualID")
linker_wgs[, projid := as.character(projid)]

filler_zeroes <- Vectorize(function(id, maxchar = 8) {
  paste(rep_len(0, maxchar - nchar(id)), collapse = "")
}, vectorize.args = "id")

linker_wgs[nchar(projid) < 8, 
           paste("Number if IDs with < 8 chars:", .N)]

linker_wgs[, projid := paste0(filler_zeroes(projid), projid)]

np <- np[!is.na(age_death) & !is.na(ceradsc)]

np_wgs <- merge(np, linker_wgs[, .(projid)], "projid")
