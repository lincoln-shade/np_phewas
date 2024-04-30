library(data.table)

dups <- fread("data/mega/duplicates_bt_cohorts.con")
remove_ids <- dups[, .(FID2, ID2)]

fwrite(remove_ids, 
       file = "data/mega/dup_ids.txt", 
       quote = FALSE,
       col.names = FALSE, 
       sep = ' ')
