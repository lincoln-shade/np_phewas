#=====================
# Make phenotype file
#=====================

library(data.table)
library(readxl)

rosmap <- as.data.table(read_excel(
  "/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))

ids <- fread("data/rosmap/rosmap_np.fam")
setnames(ids, c('V1', 'V2'), c('FID', 'IID'))
ids[, IID := as.character(IID)]

rosmap_np <- merge(ids[, .(FID, IID)], 
                   rosmap, 
                   by.x = 'IID',
                   by.y = 'projid')

ordinal_vars <- c('arteriol_scler', 'cvda_4gp2', 'tdp_st4', 'dlbdx',
                  'ceradsc', 'niareagansc', 'caa_4gp', 'braaksc')
binary_vars <- c('ci_num2_gct', 'ci_num2_mct', 'hspath_typ')
log_transform_vars <- c('nft', 'tangles', 'plaq_d', 'amyloid', 'gpath')
keep_vars <- c('FID', 'IID', binary_vars, ordinal_vars, log_transform_vars)

# dichotomize ordinal variables
# using RADC codebook recommendations if present
rosmap_np[!is.na(arteriol_scler), 
          arteriol_scler_bin := ifelse(arteriol_scler < 2, 0, 1)]

rosmap_np[!is.na(cvda_4gp2), 
          cvda_4gp2_bin := ifelse(cvda_4gp2 < 2, 0, 1)]

rosmap_np[!is.na(tdp_st4), 
          tdp_st4_bin := ifelse(tdp_st4 < 2, 0, 1)] # RADC rec

rosmap_np[!is.na(dlbdx), 
          dlbdx_bin := ifelse(dlbdx < 3, 0, 1)] # only neocort ass. w dementia

rosmap_np[!is.na(ceradsc), 
          ceradsc_bin := ifelse(ceradsc < 3, 0, 1)] # RADC rec

rosmap_np[!is.na(niareagansc), 
          niareagansc_bin := ifelse(niareagansc > 2, 0, 1)] # RADC rec

rosmap_np[!is.na(caa_4gp), 
          caa_4gp_bin := ifelse(caa_4gp < 2, 0, 1)] 

rosmap_np[!is.na(braaksc), 
          braaksc56 := ifelse(braaksc < 5, 0, 1)] 
