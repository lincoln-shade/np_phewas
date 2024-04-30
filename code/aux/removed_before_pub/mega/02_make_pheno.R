library(data.table)
library(stringi)
library(readxl)

adc <- fread("data/adc/adc_np.pheno", na.strings = '-1')
ros <- fread("data/rosmap/rosmap_np.pheno", na.strings = '-1')
act <- fread("data/act/act_np.pheno", na.strings = '-1')
adni <- fread("data/adni/adni_np.pheno", na.strings = '-1')
fam <- fread("data/mega/mega_np.fam")

act_np <- as.data.table(read_xlsx("raw_data/ACT_235_provided_2023_06_22/np.xlsx",
                                  sheet = 1))
adc_np <- setDT(readRDS("data/adc/np_qced.Rds")) # %>%
  # replace_with_na_all(condition = ~.x == -1)
setDT(adc_np)
ros_np <- setDT(read_xlsx(
  "/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))
adni_np <- fread("raw_data/NEUROPATH_07_06_21.csv")
adni_fam <- fread("data/adni/adni_np.fam")
adni_np <- merge(adni_fam, adni_np, by.x = 'V2', by.y = 'RID')
rm(adni_fam)


set_ids_as_char <- function(dt) {
  dt[, IID := as.character(IID)]
  dt[, FID := as.character(FID)]
}

set_ids_as_char(adc)
set_ids_as_char(ros)
set_ids_as_char(act)
set_ids_as_char(adni)

ros_np[, FID := '1']
act_np[, FID := '2']
setnames(ros_np, 'projid', 'IID')
setnames(act_np, 'IDfromDave', 'IID')
setnames(adni_np, c('V1', 'V2'), c('FID', 'IID'))
set_ids_as_char(adc_np)
set_ids_as_char(ros_np)
set_ids_as_char(act_np)
set_ids_as_char(adni_np)




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

#-----------------------------------------------------------------------------
# Merge all the study files, keeping variables to create additional NPE vars
#-----------------------------------------------------------------------------
mega_np <- rbindlist(list(adc_np, ros_np, act_np, adni_np), fill = TRUE)
mega_np <- merge(mega_np, pheno, c('FID', 'IID'))
setnames(mega_np, 'braak56.y', 'braak56') # braak56 also ACT variable

# braak score ordinal
mega_np[!is.na(NACCBRAA), braak := NACCBRAA] # NACC
mega_np[!is.na(braaksc), braak := as.integer(braaksc)] # ROSMAP
mega_np[!is.na(NPBRAAK), braak := NPBRAAK] # ADNI (braak existing ACT var)

# neuritic plaque CERAD score ordinal
mega_np[!is.na(NACCNEUR), cerad := NACCNEUR] # NACC
mega_np[ceradsc == 4, cerad := 0] # ROSMAP
mega_np[ceradsc == 3, cerad := 1] # ROSMAP
mega_np[ceradsc == 2, cerad := 2] # ROSMAP
mega_np[ceradsc == 1, cerad := 3] # ROSMAP
mega_np[!is.na(NPNEUR), cerad := NPNEUR] # ADNI (cerad existing ACT var)

# Thal CERAD score / diffuse plaques ordinal and binary vars
mega_np[!is.na(NACCDIFF), diffuse_abeta := NACCDIFF] # NACC/ADC
mega_np[!is.na(thal_a), diffuse_abeta := thal_a] # ACT
mega_np[!is.na(NPDIFF), diffuse_abeta := NPDIFF] # ADNI
mega_np[plaq_d == 0, diffuse_abeta := 0] # ROSMAP
mega_np[plaq_d > 0 & plaq_d <= 0.5, diffuse_abeta := 1] # ROSMAP
mega_np[plaq_d > 0.5 & plaq_d <= 1, diffuse_abeta := 2] # ROSMAP
mega_np[plaq_d > 1, diffuse_abeta := 3] # ROSMAP
mega_np[diffuse_abeta %in% 0:2, diffuse_abeta3 := 0]
mega_np[diffuse_abeta %in% 3, diffuse_abeta3 := 1]

# TDP-43 ordinal and binary vars
# Stage 1: Amygdala only
# Stage 2: + Hippocampus
# Stage 3: + Middle frontal gyrus or cortex
mega_np[NPTDPB == 0 | NPTDPC == 0 | NPTDPD == 0 | NPTDPE == 0, 
        late := 0] # NACC and ADNI
mega_np[NPTDPB == 1, late := 1] # NACC and ADNI
mega_np[NPTDPC == 1 | NPTDPD == 1, late := 2] # NACC and ADNI
mega_np[NPTDPE == 1, late := 3] # NACC and ADNI
mega_np[!is.na(tdp_st4), late := tdp_st4] # ROSMAP
mega_np[late %in% 0:1, late23 := 0]
mega_np[late %in% 2:3, late23 := 1]

# arteriolosclerosis
mega_np[!is.na(NACCARTE), arteriol := NACCARTE] # NACC
mega_np[!is.na(arteriol_scler), arteriol := arteriol_scler] # ROSMAP
mega_np[!is.na(NPARTER), arteriol := NPARTER] # ADNI
mega_np[!is.na(micro_arteriolosclerosis_id), 
        arteriol := micro_arteriolosclerosis_id - 1] # ACT

# atherosclerosis in circle of willis
mega_np[NACCAVAS %in% 0:3, atheroscler := NACCAVAS]
mega_np[cvda_4gp2 %in% 0:3, atheroscler := cvda_4gp2]
mega_np[ge_atherosclerosis_id %in% 0:3, atheroscler := ge_atherosclerosis_id]
mega_np[NPAVAS %in% 0:3, atheroscler := NPAVAS]

# lewy body pathology
mega_np[NACCLEWY %in% 0:3, lewy_body := NACCLEWY] # NACC
mega_np[NACCLEWY == 4 & NPLBOD == 5, lewy_body := 0] # NACC olfactory
mega_np[dlbdx %in% 0:3, lewy_body := dlbdx] # ROSMAP
mega_np[FID == "ADNI" & NPLBOD %in% 0:3, lewy_body := NPLBOD] # ADNI
mega_np[FID == "ADNI" & NPLBOD == 4, lewy_body := 2] # ADNI amygdala
mega_np[FID == "ADNI" & NPLBOD == 5, lewy_body := 0] # ADNI olfactory
mega_np[lbsubsnigra == 0 & lblocuscer == 0 & lbamygdala == 0 & lbfrontcor == 0,
        lewy_body := 0] # ACT
mega_np[(lbsubsnigra == 1 | lblocuscer == 1) & 
          (lbamygdala == 0 & lbfrontcor == 0), lewy_body := 1]
mega_np[(lbamygdala == 1 & lbfrontcor == 0),lewy_body := 2]
mega_np[lbfrontcor == 1, lewy_body := 3]

mega_np[lewy_body %in% 0, lewy_body123 := 0]
mega_np[lewy_body %in% 1:3, lewy_body123 := 1]

# cerebral amyloid angiopathy
mega_np[NACCAMY %in% 0:3, caa_ord := NACCAMY]
mega_np[micro_amyloidangiopathyoccipital %in% 1:4, 
        caa_ord := micro_amyloidangiopathyoccipital - 1]
mega_np[caa_4gp %in% 0:3, caa_ord := caa_4gp]
mega_np[NPAMY %in% 0:3, caa_ord := NPAMY]


plink_pheno_vars <- c("FID", "IID", "hs", "microinf", "grossinf")

fwrite(mega_np[, ..plink_pheno_vars], 
       file = "data/mega/mega_np.pheno", 
       quote = FALSE, 
       sep = ' ', 
       na = -1)

fwrite(mega_np[, .(FID, IID, braak, cerad, diffuse_abeta, arteriol, 
                   atheroscler, late, lewy_body, caa_ord)],
       file = "data/mega/mega_np_ord.pheno",
       quote = FALSE,
       sep = ' ',
       na = -1)
saveRDS(mega_np, file = "data/mega/mega_np.Rds")

# cohort-stratified pheno files
## nacc/adc
fwrite(mega_np[FID == "0", ..plink_pheno_vars], 
       file = "data/mega/nacc_np.pheno", 
       quote = FALSE, 
       sep = ' ', 
       na = -1)

fwrite(mega_np[FID == "0", 
               .(FID, IID, braak, cerad, diffuse_abeta, arteriol, 
                 atheroscler, late, lewy_body, caa_ord)],
       file = "data/mega/nacc_np_ord.pheno",
       quote = FALSE,
       sep = ' ',
       na = -1)

## rosmap
fwrite(mega_np[FID == "1", ..plink_pheno_vars], 
       file = "data/mega/rosmap_np.pheno", 
       quote = FALSE, 
       sep = ' ', 
       na = -1)

fwrite(mega_np[FID == "1", 
               .(FID, IID, braak, cerad, diffuse_abeta, arteriol, 
                 atheroscler, late, lewy_body, caa_ord)],
       file = "data/mega/rosmap_np_ord.pheno",
       quote = FALSE,
       sep = ' ',
       na = -1)

## act (no LATE data available)
fwrite(mega_np[FID == "2", ..plink_pheno_vars], 
       file = "data/mega/act_np.pheno", 
       quote = FALSE, 
       sep = ' ', 
       na = -1)

fwrite(mega_np[FID == "2", 
               .(FID, IID, braak, cerad, diffuse_abeta, arteriol, 
                 atheroscler, lewy_body, caa_ord)],
       file = "data/mega/act_np_ord.pheno",
       quote = FALSE,
       sep = ' ',
       na = -1)

## adni
fwrite(mega_np[FID == "ADNI", ..plink_pheno_vars], 
       file = "data/mega/adni_np.pheno", 
       quote = FALSE, 
       sep = ' ', 
       na = -1)

fwrite(mega_np[FID == "ADNI", 
               .(FID, IID, braak, cerad, diffuse_abeta, arteriol, 
                 atheroscler, late, lewy_body, caa_ord)],
       file = "data/mega/adni_np_ord.pheno",
       quote = FALSE,
       sep = ' ',
       na = -1)

