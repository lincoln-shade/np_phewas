library(data.table)
library(magrittr)
library(stringi)
library(readxl)

na_strings <- c("-9", "-4", "9999", "999", "888")
# Universal Data Set + NP
uds <- fread("/data_global/nacc/investigator_nacc57.csv", 
             na.strings = na_strings)
# NP data not longitudinal, so just keep obs from last visit
uds <- uds[NACCVNUM == NACCAVST]

# Minimal Data Set + NP
mds <- fread("/data_global/nacc/fardo09062019.csv", 
             na.strings = na_strings)

#------------------------------------------------------
# subset participants with NP data and NP variables
#------------------------------------------------------

# all autopsied have NPFORMVER
uds[!(is.na(NPFORMVER)), np := 1]
table(uds$np, useNA = "a")

mds[!(is.na(NPFORMVER)), np := 1]
table(mds$np, useNA = "a")

np_vars <- c(colnames(uds)[grep("NP", colnames(uds))])
np_vars_nacc <- fread("data/nacc_np_vars_names_nacc.txt", header = FALSE)$V1
np_vars_remove <- fread("data/nacc_uds_var_names_np.txt", header = FALSE)$V1
np_vars <- np_vars[!(np_vars) %in% np_vars_remove]
np_vars <- c(np_vars_nacc, np_vars)
uds_np <- uds[np == 1, c(..np_vars, "RACE", "HISPANIC"), with =FALSE]
uds_np[, uds := 1]
# month of death and year of death not in MDS, but still want to have for UDS
mds[, `:=`(NACCMOD = NA, NACCYOD = NA)]
np_mds_vars <- np_vars[!(np_vars %in% 
                           c("NPARTAG", "NPATGSEV", "NPATGAMY",
                             "NPATGAM1", "NPATGAM2", "NPATGAM3", 
                             "NPATGAM4", "NPATGAM5", "NPATGFRN", 
                             "NPATGFR1", "NPATGFR2", "NPATGFR3", 
                             "NPATGFR4"))]
mds_np <- mds[np == 1, c(np_mds_vars, "MDSRACE", "MDSHISPANIC"), 
              with = FALSE]
mds_np[, uds := 0]

np <- rbind(uds_np, mds_np, fill=TRUE)
# one duplicate ID, same values in MDS and UDS so I'll keep UDS row
sum(duplicated(np$NACCID))
np_dup_id <- np[NACCID %in% NACCID[duplicated(NACCID)], 
                -c("NACCMOD", "NACCYOD", "uds")]
np <- np[!((NACCID %in% np_dup_id$NACCID) & uds == 0)]

# variables that won't be treated as outcome variables. 
# May be used as covariates in some/all analyses
not_outcome_vars <- c("NACCID", "NACCADC", "NPFORMVER", "NPSEX", "NACCDAGE", 
                      "NACCMOD", "NACCYOD", "NACCINT", "NPPMIH", "NPFIX", 
                      "NPFIXX", "NPWBRF", "NPTAN", "NPTANX", "NPABAN", 
                      "NPABANX", "NPASAN", "NPASANX", "NPTDPAN", "NPTDPANX", 
                      "NPHISMB", "NPHISG", "NPHISSS", "NPHIST", "NPHISO", 
                      "NPHISOX", "NPPATHOX", "NACCWRI1", "NACCWRI2", 
                      "NACCWRI3", "NACCBNKF", "NPBNKB", "NACCFORM", 
                      "NACCPARA", "NACCCSFP", "NPBNKF", "NPFAUT", "NPFAUT1", 
                      "NPFAUT2", "NPFAUT3", "NPFAUT4")

# continuous variables
cont_vars <- c("NPWBRWT")

# count variables (probably wont be used as outcomes due to low sample sizes)
count_vars <- c("NPINF1A", "NPINF1B", "NPINF1D", "NPINF1F", "NPINF2A", 
                "NPINF2B", "NPINF2D", "NPINF2F", "NPINF3A", "NPINF3B", 
                "NPINF3D", "NPINF3F", "NPINF4A", "NPINF4B", "NPINF4D", 
                "NPINF4F")

# binary variables
#-------------------
# note that 0 in NACCBRNN means that NO NP change is present
bin_vars_0_1_8 <- c("NACCBRNN")
bin_vars_0_1_8_9_neg4 <- c("NPGRLA", "NPINF", "NPHEMO", "NPHEMO1", "NPHEMO2", 
                           "NPHEMO3", "NPOLD", "NPOLDD", "NPPATH", "NPPATH2", 
                           "NPPATH3", "NPPATH4", "NPPATH5", "NPPATH6", 
                           "NPPATH7", "NPPATH8", "NPPATH9", "NPPATH10", 
                           "NPPATH11", "NPFTDTAU", "NPFTDT2", "NPFTDT5", 
                           "NPFTDT6", "NPFTDT7", "NPFTDT8", "NPFTDT9", 
                           "NPFTDT10", "NPFTDTDP", "NPOFTD", "NPOFTD1", 
                           "NPOFTD2", "NPOFTD3", "NPOFTD4", "NPOFTD5", 
                           "NPTDPA", "NPTDPB", "NPTDPC", "NPTDPD", "NPTDPE", 
                           "NPPDXA", "NPPDXB", "NPPDXD", "NPPDXE", "NPPDXF", 
                           "NPPDXG", "NPPDXH", "NPPDXI", "NPPDXJ", "NPPDXK", 
                           "NPPDXL", "NPPDXM", "NPPDXN", "NPPDXP", "NPPDXQ")
bin_vars_0_1_9 <- c("NACCVASC", "NACCINF")
bin_vars_0_1_8_9 <- c("NACCMICR", "NACCHEM", "NACCNEC", "NACCPICK", "NACCCBD", 
                      "NACCPROG", "NACCPRIO", "NACCOTHP")
bin_vars_0_1_neg4 <- c("NPPATHO")
# 3 is missing
# 1 = yes, 2 = no
bin_vars_1_2_3_9_neg4 <- c("NPLINF", "NPLAC", "NPHEM", "NPMICRO", "NPART", 
                           "NPOANG", "NPSCL", "NPFRONT", "NPTAU", "NPFTDNO", 
                           "NPFTDSPC")
# 1 = yes
bin_vars_1_7 <- c("NACCDOWN")

# ordinal vars
#----------------
ord_vars_0_1_2_3_8_9_neg4 <- 
  c("NPGRCCA", "NPGRHA", "NPGRSNH", "NPGRLCH", 
    "NPADNC", "NPOLD1", "NPOLD2", "NPOLD3", "NPOLD4", "NPOLDD1", "NPOLDD2", 
    "NPOLDD3", "NPOLDD4", "NPWMR", "NPNLOSS", "NPHIPSCL")
ord_vars_0_1_2_3_8_9 <- c("NACCAVAS", "NACCNEUR", "NACCDIFF", "NACCAMY", 
                          "NACCARTE")
ord_vars_0_1_2_3_4_5_8_9_neg4 <- c("NPTHAL")
# 7 should be missing in NACCBRAA bc indicates some other tauopathy
ord_vars_0_1_2_3_4_5_6_7_8_9_neg4 <- c("NACCBRAA")

# categorical unordered variables
#---------------------------------
cat_vars_0_1_2_3_4_8_9 <- c("NACCLEWY")
cat_vars_0_1_2_3_4_5_8_9_neg4 <- c("NPLBOD", "NPALSMND")
# 3 = none present, 4 = missing
cat_vars_1_2_3_4_9_neg4 <- c("NPFTD")

#------------------------------------
# functions to set values to missing
#------------------------------------
# note, for PLINK, default missing phenotype value is -9, 
# can be explicitly noted with --missing-phenotype <integer>
# --1 flag can be used to denote that phenotypes are 0 = control, 1 = case
set_missing_vals <- function(var, na_vals) {
  var <- ifelse(var %in% na_vals, NA, var)
  return(var)
} %>% 
  Vectorize()

set_missing_8_9_neg4 <- function(var, na_vals = c(8, 9, -4)) {
  set_missing_vals(var, na_vals)
}

#-------------------------------
# set values to missing
#-------------------------------

vars_set_8_9_neg4_na <- c(bin_vars_0_1_8, 
                          bin_vars_0_1_8_9, 
                          bin_vars_0_1_8_9_neg4, 
                          bin_vars_0_1_9, 
                          bin_vars_0_1_neg4, 
                          ord_vars_0_1_2_3_8_9, 
                          ord_vars_0_1_2_3_8_9_neg4, 
                          ord_vars_0_1_2_3_4_5_8_9_neg4)

for (i in vars_set_8_9_neg4_na) {
  # np[[i]] <- set_missing_8_9_neg4(np[[i]])
  np[, (i) := set_missing_8_9_neg4(get(i))]
}

# these variables also need to have values == 2 set to 0
for (i in bin_vars_1_2_3_9_neg4) {
  np[, (i) := set_missing_vals(get(i), na_vals = c(3, 9, -4))]
  np[get(i) == 2, (i) := 0]
}

# braak variable also needs 7 set to NA
np[, NACCBRAA := set_missing_vals(NACCBRAA, c(7, 8, 9, -4))]

np <- np[!is.na(NACCNEUR)]

np[, RACE := set_missing_vals(RACE, c(50, 99))]
np[, MDSRACE := set_missing_vals(MDSRACE, c(50, 99))]
np[!is.na(MDSRACE) & is.na(RACE), RACE := MDSRACE]

np[, HISPANIC := set_missing_vals(HISPANIC, 9)]
np[, MDSHISPANIC := set_missing_vals(MDSHISPANIC, 9)]
np[MDSHISPANIC == 2, MDSHISPANIC := 0]
np[!is.na(MDSHISPANIC) & is.na(HISPANIC), HISPANIC := MDSHISPANIC]

#-----------------------------------
# Non-hispanic white participants
#-----------------------------------

nhw_ids <- fread("tmp/np_ids.tmp", header = F)
adc_fam <- fread("data/adc/adc.fam")
adc_fam[grep("[[:digit:]]_NACC[[:digit:]]{6}", V2)]
adc_fam[, V2 := stri_replace_first_regex(V2, ".*_", "")]
adc_fam[grep("[[:digit:]]_NACC[[:digit:]]{6}", V2)]
# one duplicate ID now formatted incorrectly, not an issue
adc_fam[grep("NACC[[:digit:]]{6}", V2, invert = T)]

np_nhw <- merge(np, nhw_ids[, .(V2)], by.x = "NACCID", by.y = "V2")
#-----------------------------------
# Black participants
#-----------------------------------

aa_adc1_2 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_AA/ADC1-2_AA/CleanedGenotypes/",
  "ADCR1_ADCGR2_merged_flipped_QCed_IBD0.41_AUTOSOMAL_Markers.no_ambiguous.fam"
  ))
aa_adc3 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_AA/ADC3_AA/CleanedGenotypes/",
  "adgcr3_aa_QCed+outliersRemoved_IBD0.41_AUTOSOMAL_MARKERS.no_ambiguous.fam"
  ))
aa_adc8 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_AA/ADC8_AA/CleanedGenotypes/",
  "ADC8_AA_nonexm_exome_clean.fam"
  ))
aa_adc9 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_AA/ADC9_AA/CleanedGenotypes/",
  "ADC9_newSEX_newDX_AA_clean.fam"
  ))
aa_adc10 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_AA/ADC10_AA/CleanedGenotypes/",
  "ADC10_newSEX_newDX_AA_clean.fam"
  ))
aa_adc11 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_AA/ADC11_AA/CleanedGenotypes/",
  "ADC11_aframr.clean.fam"
  ))
aa_adc12 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_AA/ADC12_AA/CleanedGenotypes/",
  "ADC12.clean.afr_amr.fam"
  ))

# reformat IDs in cohorts 1 and 2
aa_adc1_2[, V2 := 
            stri_replace_all_regex(
              V2, "[[:digit:]]{2}AD[[:digit:]]{4}", "")]
aa_adc1_2[, V2 := stri_replace_all_regex(V2, "_", "")]

aa_adc <- rbindlist(list(aa_adc1_2, aa_adc3, aa_adc8, aa_adc9, aa_adc10, 
                    aa_adc11, aa_adc12))

np_aa <- np[NACCID %in% c(aa_adc1_2$V2, aa_adc8$V2, aa_adc9$V2, aa_adc10$V2,
                          aa_adc11$V2, aa_adc12$V2)]

merge(merge(adc_fam, aa_adc[, .(V2)], by = "V2"), 
      np[, .(NACCID)], by.x = "V2", by.y = "NACCID")

#-------------------------------------
# Hispanic participants
#-------------------------------------
hisp_adc8 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_Hispanic/ADC8_Hispanic/CleanedGenotypes/",
  "ADC8_Hispanic_nonexm_exome_clean.fam"
))

hisp_adc9 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_Hispanic/ADC9_Hispanic/CleanedGenotypes/",
  "ADC9_newSEX_newDX_hispanic_clean.fam"
))

hisp_adc10 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_Hispanic/ADC10_Hispanic/CleanedGenotypes/",
  "ADC10_newSEX_newDX_hispanic_clean.fam"
))

hisp_adc11 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_Hispanic/ADC11_Hispanic/CleanedGenotypes/",
  "ADC11_hispanic.clean.fam"
))

hisp_adc12 <- fread(paste0(
  "/data_global/ADGC_GWAS/ADGC_Hispanic/ADC12_Hispanic/CleanedGenotypes/",
  "ADC12.clean.hispanic.fam"
))

np_hisp <- np[NACCID %in% c(hisp_adc8$V2, hisp_adc9$V2, hisp_adc10$V2,
                            hisp_adc11$V2, hisp_adc12$V2)]
# np_hisp[NACCID %in% nacc_hisp$sample, .(NACCID, NACCDAGE)]
# np[NACCID %in% c(np_hisp$NACCID, nacc_hisp$sample), .(NACCID, NACCDAGE, NPTDPE)]

#-------------------------------------------------
# ADSP participants
#-------------------------------------------------
adsp <- setDT(read_xlsx(
  "/data_global/adsp/ADSPID_to_NACCID_mapping_WES-WGS.xlsx"))

np_wgs <- np[NACCID %in% adsp[WGS == 1, NACCID], 
             .(NACCID, NACCDAGE, NACCNEUR, NPTDPE)]
np_wes <- np[NACCID %in% adsp[WES == 1, NACCID], 
             .(NACCID, NACCDAGE, NACCNEUR, NPTDPE)]
merge(np_wgs, np_wes)

np_adsp <- adsp[NACCID %in% np$NACCID, .(NACCID, WES, WGS)]

#--------------------------------------------
# GWAS, WES, and WGS data for NACC NP
#--------------------------------------------
np_gwas <- rbindlist(list(np_nhw, np_aa, np_hisp))
np_gwas[, gwas := 1]
np_gwas <- np_gwas[!duplicated(NACCID), .(NACCID, gwas)]

np_gen <- merge(np_gwas, np_adsp, all.x = T, all.y = T)
np_gen[is.na(gwas), gwas := 0]
np_gen[is.na(WES), WES := 0]
np_gen[is.na(WGS), WGS := 0]

np_gen[gwas == 1, table(WES, WGS)]
np_gen[gwas == 0, table(WES, WGS)]

np[NPTDPB == 0 | NPTDPC ==0 | NPTDPD == 0 | NPTDPE == 0, LATE := 0]
np[NPTDPC == 1 | NPTDPD == 1 | NPTDPE == 1, LATE := 1]

np_gen_late <- merge(np_gen, np[, .(NACCID, LATE)])
