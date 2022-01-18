#============================
# Make covariate file
#============================

library(data.table)
library(readxl)

rosmap <- as.data.table(read_excel(
  "/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))

load('data/rosmap/mypcair.Rdata')

ids <- fread("data/rosmap/rosmap_np.fam")
setnames(ids, c('V1', 'V2'), c('FID', 'IID'))
ids[, IID := as.character(IID)]

rosmap_np <- merge(ids[, .(FID, IID)], 
                   rosmap, 
                   by.x = 'IID',
                   by.y = 'projid')

pcs <- as.data.table(mypcair$vectors[, 1:5], keep.rownames = TRUE)

rosmap_np <- rosmap_np[, .(FID, IID, age_death, msex, study)]
rosmap_np <- merge(rosmap_np, pcs, by.x = 'IID', by.y = 'rn')
setcolorder(rosmap_np, 'FID')
setnames(rosmap_np, 'study', 'ROS')
rosmap_np[ROS == 'ROS', ROS := '1']
rosmap_np[ROS == 'MAP', ROS := '0']
rosmap_np[, table(ROS)]

fwrite(rosmap_np, 
       file = 'data/rosmap/rosmap_np.covar', 
       quote = FALSE, 
       sep = ' ')
