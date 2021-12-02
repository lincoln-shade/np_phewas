
#=================================
# create ACT .covar file
#=================================
#
# sex, age of death, ACT cohort, PCs 1-5

library(data.table)
library(readxl)
library(magrittr)

act <- as.data.table(read_xlsx("raw_data/act/DataForE235_20210421.xlsx",
                               sheet = 1))
act_np <- as.data.table(read_xlsx("raw_data/act/DataForE235_20210421.xlsx",
                                  sheet = 2))

act_cohorts <- fread("data/act_ids_adgc.txt", 
                     colClasses = c("character", "character", "factor"))
pcs <- fread("data/plink/act_np_pcs.txt")

# from https://www.maelstrom-research.org/variable/act_recruitment:gender:Collected
# gender coded as 2 female, 1 male in ACT

act_np <- merge(act_np[, .(IDfromDave, autopsyage)],
                act[, .(IDfromDave, gender)], 
                by = "IDfromDave")

act_np[gender == 2, gender := 0]
setnames(act_np, colnames(act_np), c("IID", "age_death", "male"))

act_cohort_dummy <- model.matrix(~ V3, act_cohorts)
act_cohorts <- cbind(act_cohorts, act_cohort_dummy)
act_cohorts <- act_cohorts[, .(V1, V2, V32, V33)]
setnames(act_cohorts, colnames(act_cohorts), 
         c("FID", "IID", "ACT2", "ACT3"))

act_np <- merge(act_np, 
                act_cohorts,
                by = "IID")

act_np <- merge(act_np, 
                pcs, 
                by = "IID")

act_np <- act_np[complete.cases(act_np)]
setcolorder(act_np, "FID")
act_np[, FID := 0]

act_np <- act_np[!(duplicated(IID))]

fwrite(act_np, 
       file = "data/plink/act_np.covar",
       sep = " ", 
       quote = FALSE)
