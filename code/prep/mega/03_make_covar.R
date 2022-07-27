
library(data.table)
library(stringi)

adc <- fread("data/adc/adc_np.covar")
ros <- fread("data/rosmap/rosmap_np.covar")
act <- fread("data/act/act_np.covar")
adni <- fread("data/adni/adni_np.covar")
pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = '-1')
load("data/mega/mypcair.Rdata")

ros[, rosmap := 1]
act[, act := 1]
adni[, adni := 1]

adc_covars_uniq <- c('adgc_2', 'adgc_3', 'adgc_4', 'adgc_5', 'adgc_6',
                     'adgc_7', 'adgc_8', 'adgc_9', 'adgc_10', 'adgc_11', 
                     'adgc_12')
ros_covars_uniq <- c('rosmap', 'ROS')
act_covars_uniq <- c('act', 'ACT2', 'ACT3')
adni_covars_uniq <- c('adni')

add_covars_set_zero <- function(dt, vars) {
  for (i in vars) {
    dt[, (i) := 0]
  }
}

add_covars_set_zero(adc, c(ros_covars_uniq, act_covars_uniq, adni_covars_uniq))
add_covars_set_zero(ros, c(adc_covars_uniq, act_covars_uniq, adni_covars_uniq))
add_covars_set_zero(act, c(ros_covars_uniq, adc_covars_uniq, adni_covars_uniq))
add_covars_set_zero(adni, c(adc_covars_uniq, act_covars_uniq, ros_covars_uniq))

set_pcs_null <- function(dt) {
  pc_vars <- colnames(dt)[grep('pc', colnames(dt))]
  for (pc in pc_vars) {
    dt[, (pc) := NULL]
  }
}

set_pcs_null(adc)
set_pcs_null(ros)
set_pcs_null(act)
set_pcs_null(adni)

set_ids_as_char <- function(dt) {
  dt[, IID := as.character(IID)]
  dt[, FID := as.character(FID)]
}

set_ids_as_char(adc)
set_ids_as_char(ros)
set_ids_as_char(act)
set_ids_as_char(adni)

covar <- rbindlist(list(adc, ros, act, adni), use.names = TRUE)

n_pcs <- 10
pcs <- as.data.table(mypcair$vectors[, 1:(min(32, n_pcs))], 
                     keep.rownames = TRUE)
colnames(pcs) <- stri_replace_all_fixed(colnames(pcs), 'V', 'pc')
setnames(pcs, 'rn', 'IID')

covar <- merge(covar, pcs, 'IID')
covar <- merge(covar, pheno[, .(FID, IID)], c('FID', 'IID'))

for (i in 1:n_pcs) {
  covar[, (paste0('pc', i)) := scale(get(paste0('pc', i)))]
}

covar[, age_death2 := NULL]

fwrite(covar, 
       file = "data/mega/mega_np.covar", 
       quote = FALSE, 
       sep = ' ', 
       na = NA)

# subset age 80+ at death
# remove ADNI since only 26 ADNI 80+, could affect models
fwrite(covar[age_death >= 80 & FID != "ADNI"], 
       file = "data/mega/mega_np_age80.covar", 
       quote = FALSE, 
       sep = ' ', 
       na = NA)

covar <- merge(covar, pheno[, .(FID, IID, cerad, diffuse_abeta)])

fwrite(covar, 
       file = "data/mega/mega_np_part_pos.covar", 
       quote = FALSE, 
       sep = ' ', 
       na = NA)


