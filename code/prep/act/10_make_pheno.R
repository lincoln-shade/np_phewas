#==================================
# categorize the act NP variables
# by outcome type, coding, etc.
#==================================

library(data.table)
library(readxl)
library(magrittr)

act <- setDT(read_xlsx("raw_data/ACT_235_2023_08_21/np_2023_08_21.xlsx"))

act_fam <- fread("data/act/act_np.fam") %>% 
  setnames(., c("V1", "V2"), c("FID", "IID"))

setnames(act, "IDFromDave", "IID")

act <- merge(act_fam[, .(FID, IID)], act, "IID")
setcolorder(act, "FID")

# IIDs missing cerad23 are also missing all other phenotypes, so remove those
# update as of 2023-08-23 -- no longer true, both samples missing have LATE data
# act <- act[!(is.na(cerad_score))]
# 

# rename phenotypes
setnames(
  act, 
  c("late_stg", "cerad_score", "braak_stg", "amylangi_score", "athero_level",
    "arteriolo_level", "thal_phase", "hipscl_any", "chronic_grossinfarcts_any", "chronic_microinfarcts_any"), 
  c("late", "cerad", "braak", "caa", "atheroscler", "arteriol", "diffuse_abeta", "hs", "grossinf", "microinf")
  )

# create ordinal lewy body variable from categorical variable
act[
  lewybodydis_cat == 1 | lewybodydis_cat == 1.5 | lewybodydis_cat == 5,
  lewy_body := 0
]
act[lewybodydis_cat == 3, lewy_body := 1]
act[lewybodydis_cat == 2 | lewybodydis_cat == 4, lewy_body := 2]
act[lewybodydis_cat == 6, lewy_body := 3]

keep_vars = c(
  "FID", "IID", "late", "cerad", "braak", "caa", "atheroscler", "arteriol", 
  "diffuse_abeta", "lewy_body", "hs", "grossinf", "microinf"
)

act[diffuse_abeta == 5, diffuse_abeta := 4]
act[diffuse_abeta > 0, diffuse_abeta := diffuse_abeta - 1]

fwrite(act[, ..keep_vars],
       file = "data/act/act_np.pheno",
       quote = FALSE, 
       sep = " ", 
       na = -1)
