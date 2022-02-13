library(data.table)
library(stringi)

adc <- fread("data/adc/adc_np.pheno")
ros <- fread("data/rosmap/rosmap_np.pheno")
act <- fread("data/act/act_np.pheno")
adni <- fread("data/adni/adni_np.pheno")
fam <- fread("data/mega/mega_np.fam")

set_ids_as_char <- function(dt) {
  dt[, IID := as.character(IID)]
  dt[, FID := as.character(FID)]
}

set_ids_as_char(adc)
set_ids_as_char(ros)
set_ids_as_char(act)
set_ids_as_char(adni)

# rename ADC, ROSMAP, and ACT phenotypes to match those in ADNI
# (since I renamed those myself already)
phenotypes <- colnames(adni)[3:length(colnames(adni))]
setnames(adc, c('LEWY', 'ATHCW', 'NACCINF', 'NACCMICR', 'HS', 'BRAAK', 
                'NEUR', 'CAA', 'ASC'), 
         phenotypes)
setnames(ros, 
         c('dlbdx_bin', 'cvda_4gp2_bin', 'ci_num2_gct', 'ci_num2_mct',
           'hspath_typ', 'braaksc_bin', 'ceradsc_bin', 'caa_4gp_bin',
           'arteriol_scler_bin'),
         phenotypes)
setnames(act, 
         c('any_lb4', 'ath_bin', 'any_macro', 'any_mvl', 'any_hs', 'braak56',
           'cerad3', 'caa_bin', 'asc_bin'),
         phenotypes)

keep_vars <- c('FID', 'IID', phenotypes)
pheno <- rbindlist(list(adc[, ..keep_vars], 
                        ros[, ..keep_vars], 
                        act[, ..keep_vars],
                        adni[, ..keep_vars]))

pheno <- merge(pheno, 
               fam[, .(V1, V2)], 
               by.x = c('FID', 'IID'), 
               by.y = c('V1', 'V2'))

fwrite(pheno, 
       file = "data/mega/mega_np.pheno", 
       quote = FALSE, 
       sep = ' ', 
       na = -1)

saveRDS(pheno, file = "data/mega/mega_np.Rds")
