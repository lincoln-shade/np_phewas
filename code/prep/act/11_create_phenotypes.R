#==================================
# categorize the NACC NP variables
# by outcome type, coding, etc.
#==================================

library(data.table)
library(readxl)
library(magrittr)

act <- as.data.table(read_xlsx("raw_data/act/DataForE235_20210421.xlsx",
                               sheet = 2))

act_fam <- fread("data/plink/act_np.fam") %>% 
  setnames(., c("V1", "V2"), c("FID", "IID"))

adc <- fread("data/plink/adc_np.pheno")

act_phenos_in_nacc <- c("any_macro", 
                        "calc_haas_cerebralmicroinfarcts",
                        "calc_haas_deepgraymicroinfarcts",
                        # "ap_freshbrainweight",
                        "micro_arteriolosclerosis_id",
                        "any_hs",
                        "braak56",
                        "ge_atherosclerosis_id", 
                        "cerad23")

nacc_phenos_in_act <- c("NPINF",
                        "NPOLD1",
                        "NPOLD3",
                        # "NPWBRWT", # continuous, not in pheno file
                        "ASC", 
                        "HS",
                        "BRAAK",
                        "ATHCW",
                        "NEUR")

act_pheno <- act[, c("IDfromDave", ..act_phenos_in_nacc)]

act_pheno[calc_haas_cerebralmicroinfarcts == 0, 
          calc_haas_cerebralmicroinfarcts_bin := 0]

act_pheno[calc_haas_cerebralmicroinfarcts > 0, 
          calc_haas_cerebralmicroinfarcts_bin := 1]

act_pheno[calc_haas_deepgraymicroinfarcts == 0, 
          calc_haas_deepgraymicroinfarcts_bin := 0]

act_pheno[calc_haas_deepgraymicroinfarcts > 0, 
          calc_haas_deepgraymicroinfarcts_bin := 1]

act_pheno[micro_arteriolosclerosis_id < 3, 
          micro_arteriolosclerosis_id_bin := 0]

act_pheno[micro_arteriolosclerosis_id >= 3, 
          micro_arteriolosclerosis_id_bin := 1]

act_pheno[ge_atherosclerosis_id < 3, 
          ge_atherosclerosis_id_bin := 0]

act_pheno[ge_atherosclerosis_id >= 3,
          ge_atherosclerosis_id_bin := 1]

setnames(act_pheno, "IDfromDave", "IID")

act_pheno <- act_pheno[, .(IID, 
                           any_macro, 
                           calc_haas_cerebralmicroinfarcts_bin, 
                           calc_haas_deepgraymicroinfarcts_bin, 
                           micro_arteriolosclerosis_id_bin,
                           any_hs, 
                           braak56, 
                           ge_atherosclerosis_id_bin, 
                           cerad23)]

act_pheno <- merge(act_fam[, .(FID, IID)], act_pheno, "IID")
setcolorder(act_pheno, "FID")

# IIDs missing cerad23 are also missing all other phenotypes, so remove those
act_pheno <- act_pheno[!(is.na(cerad23))]

fwrite(act_pheno,
       file = "data/plink/act_np.pheno",
       quote = FALSE, 
       sep = " ", 
       na = -1)
