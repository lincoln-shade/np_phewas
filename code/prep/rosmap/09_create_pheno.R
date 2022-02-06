#=====================
# Make phenotype file
#=====================

library(data.table)
library(magrittr)
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
# log_transform_vars <- c('nft', 'tangles', 'plaq_d', 'amyloid', 'gpath')


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
          braaksc_bin := ifelse(braaksc < 5, 0, 1)] 

# log_var <- function(x) {
#   rosmap_np[, (paste0(x, '_log')) := log(get(x) + 0.01)]
# } 
# 
# for (v in log_transform_vars) {
#   log_var(v)
# }
bin_ord_vars <- paste0(ordinal_vars, "_bin")
keep_vars <- c('FID', 'IID', binary_vars, bin_ord_vars)
rosmap_np <- rosmap_np[, ..keep_vars]
index_keep <- rosmap_np[, -c('FID', 'IID')][, rowSums(is.na(.SD)) != ncol(.SD)]
rosmap_np <- rosmap_np[index_keep]

fwrite(rosmap_np, 
       file = 'data/rosmap/rosmap_np.pheno',
       quote = FALSE,
       sep = ' ',
       na = '-1')
