#==================================
# categorize the NACC NP variables
# by outcome type, coding, etc.
#==================================

library(pacman)
p_load(data.table, magrittr)

np <- readRDS("data/nacc_np.Rds")

# variables that won't be treated as outcome variables. May be used as covariates in some/all analyses
not_outcome_vars <- c("NACCID", "NACCADC", "NPFORMVER", "NPSEX", "NACCDAGE", "NACCMOD", "NACCYOD", "NACCINT", "NPPMIH", "NPFIX", "NPFIXX", "NPWBRF", "NPTAN", 
                      "NPTANX", "NPABAN", "NPABANX", "NPASAN", "NPASANX", "NPTDPAN", "NPTDPANX", "NPHISMB", "NPHISG", "NPHISSS", "NPHIST", "NPHISO", "NPHISOX", 
                      "NPPATHOX", "NACCWRI1", "NACCWRI2", "NACCWRI3", "NACCBNKF", "NPBNKB", "NACCFORM", "NACCPARA", "NACCCSFP", "NPBNKF", "NPFAUT", "NPFAUT1", 
                      "NPFAUT2", "NPFAUT3", "NPFAUT4")

# continuous variables
cont_vars <- c("NPWBRWT")

# count variables
count_vars <- c("NPINF1A", "NPINF1B", "NPINF1D", "NPINF1F", "NPINF2A", "NPINF2B", "NPINF2D", "NPINF2F", "NPINF3A", "NPINF3B", "NPINF3D", "NPINF3F", "NPINF4A", 
                "NPINF4B", "NPINF4D", "NPINF4F")

# binary variables
#-------------------
# note that 0 in NACCBRNN means that NO NP change is present
bin_vars_0_1_8 <- c("NACCBRNN")
bin_vars_0_1_8_9_neg4 <- c("NPGRLA", "NPINF", "NPHEMO", "NPHEMO1", "NPHEMO2", "NPHEMO3", "NPOLD", "NPOLDD", "NPPATH", "NPPATH2", "NPPATH3", "NPPATH4", "NPPATH5", 
                           "NPPATH6", "NPPATH7", "NPPATH8", "NPPATH9", "NPPATH10", "NPPATH11", "NPFTDTAU", "NPFTDT2", "NPFTDT5", "NPFTDT6", "NPFTDT7", 
                           "NPFTDT8", "NPFTDT9", "NPFTDT10", "NPFTDTDP", "NPOFTD", "NPOFTD1", "NPOFTD2", "NPOFTD3", "NPOFTD4", "NPOFTD5", "NPTDPA", "NPTDPB", 
                           "NPTDPC", "NPTDPD", "NPTDPE", "NPPDXA", "NPPDXB", "NPPDXD", "NPPDXE", "NPPDXF", "NPPDXG", "NPPDXH", "NPPDXI", "NPPDXJ", "NPPDXK", 
                           "NPPDXL", "NPPDXM", "NPPDXN", "NPPDXP", "NPPDXQ")
bin_vars_0_1_9 <- c("NACCVASC", "NACCINF")
bin_vars_0_1_8_9 <- c("NACCMICR", "NACCHEM", "NACCNEC", "NACCPICK", "NACCCBD", "NACCPROG", "NACCPRIO", "NACCOTHP")
bin_vars_0_1_neg4 <- c("NPPATHO")
# 3 is missing
# 1 = yes, 2 = no
bin_vars_1_2_3_9_neg4 <- c("NPLINF", "NPLAC", "NPHEM", "NPMICRO", "NPART", "NPOANG", "NPSCL", "NPFRONT", "NPTAU", "NPFTDNO", "NPFTDSPC")
# 1 = yes
bin_vars_1_7 <- c("NACCDOWN")

# ordinal vars
#----------------
ordinal_vars_0_1_2_3_8_9_neg4 <- c("NPGRCCA", "NPGRHA", "NPGRSNH", "NPGRLCH", "NPADNC", "NPOLD1", "NPOLD2", "NPOLD3", "NPOLD4", "NPOLDD1", "NPOLDD2", 
                                   "NPOLDD3", "NPOLDD4", "NPWMR", "NPNLOSS", "NPHIPSCL")
ordinal_vars_0_1_2_3_8_9 <- c("NACCAVAS", "NACCNEUR", "NACCDIFF", "NACCAMY", "NACCARTE")
ord_vars_0_1_2_3_4_5_8_9_neg4 <- c("NPTHAL")
# 7 should be missing in NACCBRAA bc indicates some other tauopathy
ord_vars_0_1_2_3_4_5_6_7_8_9_neg4 <- c("NACCBRAA")

# categorical unordered variables
#---------------------------------
cat_vars_0_1_2_3_4_8_9 <- c("NACCLEWY")
cat_vars_0_1_2_3_4_5_8_9_neg4 <- c("NPLBOD", "NPALSMND")
# 3 = none present, 4 = missing
cat_vars_1_2_3_4_9_neg4 <- c("NPFTD")