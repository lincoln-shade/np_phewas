
#------------------------------------------------
# list of ROSMAP IDs with age of death available
#------------------------------------------------

library(data.table)
library(readxl)

rosmap <- as.data.table(read_excel(
    "/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))

rosmap <- rosmap[!(is.na(age_death))]

fwrite(rosmap[, .(projid, projid)], 
       file = "data/rosmap/np_ids.txt",
       col.names = FALSE,
       quote = FALSE,
       sep = ' ')
