#============================
# Make covariate file
#============================

library(data.table)
library(readxl)
library(stringi)

rosmap <- as.data.table(read_excel(
  "/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))

load('data/rosmap/mypcair.Rdata')
n_pcs <- 10
pcs <- as.data.table(mypcair$vectors[, 1:(min(32, n_pcs))], 
                     keep.rownames = TRUE)
colnames(pcs) <- stri_replace_all_fixed(colnames(pcs), 'V', 'pc')
setnames(pcs, 'rn', 'IID')

ids <- fread("data/rosmap/rosmap_np.fam")
setnames(ids, c('V1', 'V2'), c('FID', 'IID'))
ids[, IID := as.character(IID)]

rosmap_np <- merge(ids[, .(FID, IID)], 
                   rosmap, 
                   by.x = 'IID',
                   by.y = 'projid')

rosmap_np <- rosmap_np[, .(FID, IID, age_death, msex, study)]
rosmap_np <- merge(rosmap_np, pcs, by= 'IID')
setcolorder(rosmap_np, 'FID')
setnames(rosmap_np, 'study', 'ROS')
rosmap_np[ROS == 'ROS', ROS := '1']
rosmap_np[ROS == 'MAP', ROS := '0']
rosmap_np[, table(ROS)]
rosmap_np[, age_death2 := age_death^2]


fwrite(rosmap_np, 
       file = 'data/rosmap/rosmap_np.covar', 
       quote = FALSE, 
       sep = ' ',
       na = -1)
