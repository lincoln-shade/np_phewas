
library(data.table)
library(readxl)
act_np <- as.data.table(read_xlsx("raw_data/act/DataForE235_20210421.xlsx", sheet = 2))
act_np[, FID := 0]
fwrite(act_np[, .(FID, IDfromDave)], 
       file = "data/act_np_ids.txt", 
       quote = FALSE, 
       col.names = FALSE, 
       sep = " ")
