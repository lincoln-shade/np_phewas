#==================================
# categorize the NACC NP variables
# by outcome type, coding, etc.
#==================================

library(pacman)
p_load(data.table, magrittr)

np <- readRDS("data/np.Rds")
covar <- fread("data/plink/adc_np.covar")
np <- np[IID %in% covar$IID]

apply_exclusion_criteria <- function(dt) {
  dt %>% 
    .[is.na(NACCDOWN) | NACCDOWN != 1] %>%
    .[is.na(NPPDXA) | NPPDXA != 1] %>% .[is.na(NPPDXB) | NPPDXB != 1] %>%
    .[is.na(NPPDXD) | NPPDXD != 1] %>% .[is.na(NPPDXE) | NPPDXE != 1] %>%
    .[is.na(NPPDXF) | NPPDXF != 1] %>% .[is.na(NPPDXG) | NPPDXG != 1] %>%
    .[is.na(NPPDXH) | NPPDXH != 1] %>% .[is.na(NPPDXI) | NPPDXI != 1] %>%
    .[is.na(NPPDXJ) | NPPDXJ != 1] %>% .[is.na(NPPDXK) | NPPDXK != 1] %>% 
    .[is.na(NPPDXL) | NPPDXL != 1] %>% .[is.na(NPPDXM) | NPPDXM != 1] %>% 
    .[is.na(NPPDXN) | NPPDXN != 1] %>% .[is.na(NACCPRIO) | NACCPRIO != 1] %>% 
    .[is.na(NPPATH10) | NPPATH10 != 1] %>% 
    .[is.na(NPALSMND) | NPALSMND != 1] %>% 
    .[is.na(NPFTDTAU) | NPFTDTAU != 1] %>% 
    .[is.na(NPOFTD) | NPOFTD != 1] %>% 
    .[is.na(NPFTDTDP) | NPFTDTDP != 1] %>% 
    as.data.table()
}

exclusion_vars <- c("NACCDOWN", "NPPDXA", "NPPDXB", "NPPDXD", "NPPDXE", 
                    "NPPDXF", "NPPDXG", "NPPDXH", "NPPDXI", "NPPDXJ", 
                    "NPPDXK", "NPPDXL", "NPPDXM", "NPPDXN", "NACCPRIO", 
                    "NPPATH10", "NPALSMND", "NPFTDTAU", "NPOFTD", "NPFTDTDP")

np <- apply_exclusion_criteria(np)
#----------------------------
# categorize NP variables
#----------------------------

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

#---------------------
# complex/derived NPE
#---------------------

quick_table <- function(phenotype) {
  np[, table(get(phenotype), useNA = "a")]
}
# HS binary
np[NPHIPSCL == 0 | NPSCL == 0, HS := 0]
np[NPHIPSCL %in% 1:3 | NPSCL == 1, HS := 1]
# np[NPHIPSCL == 2, HS := 2]

# B-ASC binary
np[NACCARTE %in% 0:1, ASC := 0]
np[NACCARTE %in% 2:3, ASC := 1]

# LATE
np[NPTDPB == 0 | NPTDPC ==0 | NPTDPD == 0 | NPTDPE == 0, LATE := 0]
np[NPTDPB == 1 | NPTDPC == 1 | NPTDPD == 1 | NPTDPE == 1, LATE := 1]

# Atherosclerosis in Circle of Willis
np[NACCAVAS %in% 0, ATHCW := 0]
np[NACCAVAS %in% 1:3, ATHCW := 1]

# Braak stage (B score)
np[NACCBRAA %in% 0:4, BRAAK := 0]
np[NACCBRAA %in% 5:6, BRAAK := 1]

# Density of neocorical neuritic plaques (C score)
np[NACCNEUR %in% 0:2, NEUR := 0]
np[NACCNEUR %in% 3, NEUR := 1]

# Density of diffuse plaques
np[NACCDIFF %in% 0:2, DIFF := 0]
np[NACCDIFF %in% 3, DIFF := 1]

# CAA
np[NACCAMY %in% 0, CAA := 0]
np[NACCAMY %in% 1:3, CAA := 1]

# Lewy Bodies
np[NACCLEWY %in% 0, LEWY := 0]
np[NACCLEWY %in% 1:4, LEWY := 1]

# dichotomizing various ordinal variables
ord_vars_0_1_2_3_8_9_neg4_bin <- paste0(ord_vars_0_1_2_3_8_9_neg4, "_bin")

for (i in ord_vars_0_1_2_3_8_9_neg4) {
  np[ get(i) %in% 0:1, (paste0(i, "_bin")) := 0]
  np[ get(i) %in% 2:3, (paste0(i, "_bin")) := 1]
}

#-------------------
# create new object
#-------------------

# qced data set
saveRDS(np, file = "data/np_qced.Rds")

# binary outcome variables for analysis
pheno_vars <- c("FID", "IID", "ASC", "HS", "LATE", "ATHCW", "BRAAK", 
                "DIFF", "CAA", "LEWY", 
                bin_vars_0_1_8_9_neg4, 
                bin_vars_0_1_8_9, 
                bin_vars_0_1_9, 
                bin_vars_0_1_neg4, 
                bin_vars_1_2_3_9_neg4, 
                ord_vars_0_1_2_3_8_9_neg4_bin)

fwrite(np[, ..pheno_vars][, -..exclusion_vars], 
       file = "data/plink/adc_np.pheno",
       quote = FALSE,
       sep = " ",
       na = -1)

# ordinal outcome variables for analysis
ordinal_vars <- c("FID", "IID", 
                  ord_vars_0_1_2_3_8_9_neg4, 
                  ord_vars_0_1_2_3_8_9, 
                  ord_vars_0_1_2_3_4_5_8_9_neg4, 
                  ord_vars_0_1_2_3_4_5_6_7_8_9_neg4
)

fwrite(np[, ..ordinal_vars][, -..exclusion_vars],
       file = "data/adc_np_ord.txt",
       quote = FALSE, 
       sep = " ", 
       na = -1)