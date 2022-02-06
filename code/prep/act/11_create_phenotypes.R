#==================================
# categorize the NACC NP variables
# by outcome type, coding, etc.
#==================================

library(data.table)
library(readxl)
library(magrittr)

act <- as.data.table(read_xlsx("raw_data/act/DataForE235_20210421.xlsx",
                               sheet = 2))

act_fam <- fread("data/act/act_np.fam") %>% 
  setnames(., c("V1", "V2"), c("FID", "IID"))

adc <- fread("data/adc/adc_np.pheno")

# acts_in_nacc <- c("any_macro", 
#                         "calc_haas_cerebralmicroinfarcts",
#                         "calc_haas_deepgraymicroinfarcts",
#                         "micro_arteriolosclerosis_id",
#                         "any_hs",
#                         "braak56",
#                         "ge_atherosclerosis_id", 
#                         "cerad23", 
#                         'any_lb4')
# 
# nacc_phenos_in_act <- c("NPINF",
#                         "NPOLD1",
#                         "NPOLD3",
#                         "ASC", 
#                         "HS",
#                         "BRAAK",
#                         "ATHCW",
#                         "NEUR",
#                         'NACCLEWY')

# brain arterioloscerosis binary None or Mild vs. Moderate or Severe
act[micro_arteriolosclerosis_id < 3, 
          asc_bin := 0]

act[micro_arteriolosclerosis_id >= 3, 
          asc_bin := 1]

# Atherosclerosis binary None or Mild vs. Moderate or Severe
act[ge_atherosclerosis_id < 3, 
          ath_bin := 0]

act[ge_atherosclerosis_id >= 3,
          ath_bin := 1]

# Cerebral amyloid angiopathy binary None vs. Mild, Moderate, or Severe
act[micro_amyloidangiopathyoccipital == 1, caa_bin := 0]
act[micro_amyloidangiopathyoccipital %in% 2:4, caa_bin := 1]

act <- act[, .(IDfromDave, any_hs, any_mvl, any_macro, any_lb4, asc_bin, 
               ath_bin, cerad23, braak56, caa_bin)]
setnames(act, "IDfromDave", "IID")

act <- merge(act_fam[, .(FID, IID)], act, "IID")
setcolorder(act, "FID")

# IIDs missing cerad23 are also missing all other phenotypes, so remove those
act <- act[!(is.na(cerad23))]

fwrite(act,
       file = "data/act/act_np.pheno",
       quote = FALSE, 
       sep = " ", 
       na = -1)
