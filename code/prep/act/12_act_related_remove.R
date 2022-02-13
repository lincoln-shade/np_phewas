library(data.table)

related <- fread("data/act/act_rosmap_adc.kin0")
remove_ids <- related[FID2 == 2, .(FID2, ID2)]

fwrite(remove_ids, 
       file = "data/act/act_related.remove", 
       quote = FALSE,
       col.names = FALSE, 
       sep = ' ')
