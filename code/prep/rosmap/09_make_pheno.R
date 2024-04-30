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

ordinal_vars <- c('arteriol_scler', 'cvda_4gp2', 'tdp_st4', 'dlbdx', 'caa_4gp', 'braaksc')
ordinal_vars_new <- c("arteriol", 'atheroscler', 'late', 'lewy_body', 'caa', 'braak')

binary_vars <- c('ci_num2_gct', 'ci_num2_mct', 'hspath_typ')
binary_vars_new <- c('grossinf', 'microinf', 'hs')

setnames(rosmap_np, ordinal_vars, ordinal_vars_new)
setnames(rosmap_np, binary_vars, binary_vars_new)

# cerad score low-to-high
rosmap_np[ceradsc == 4, cerad := 0]
rosmap_np[ceradsc == 3, cerad := 1]
rosmap_np[ceradsc == 2, cerad := 2]
rosmap_np[ceradsc == 1, cerad := 3]
rosmap_np[, ceradsc := NULL]

# diffuse plaques score
rosmap_np[plaq_d == 0, diffuse_abeta := 0]
rosmap_np[plaq_d > 0 & plaq_d <= 0.5, diffuse_abeta := 1]
rosmap_np[plaq_d > 0.5 & plaq_d <= 1, diffuse_abeta := 2]
rosmap_np[plaq_d > 1, diffuse_abeta := 3]
rosmap_np[, plaq_d := NULL]

ordinal_vars_new = c(ordinal_vars_new, "cerad", "diffuse_abeta")
keep_vars <- c('FID', 'IID', binary_vars_new, ordinal_vars_new)
rosmap_np <- rosmap_np[, ..keep_vars]

index_keep <- rosmap_np[, -c('FID', 'IID')][, rowSums(is.na(.SD)) != ncol(.SD)]
rosmap_np <- rosmap_np[index_keep]

fwrite(rosmap_np, 
       file = 'data/rosmap/rosmap_np.pheno',
       quote = FALSE,
       sep = ' ',
       na = '-1')
