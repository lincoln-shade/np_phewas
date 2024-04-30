
library(data.table)
library(readxl)
act_np <- setDT(read_xlsx("raw_data/ACT_235_2023_08_21/np_2023_08_21.xlsx"))
act_np[, FID := 0]
fwrite(act_np[, .(FID, IDFromDave)], 
       file = "data/act_np_ids.txt", 
       quote = FALSE, 
       col.names = FALSE, 
       sep = " ")
