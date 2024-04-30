
#=================================
# create ACT .covar file
#=================================

# sex, age of death, ACT cohort, PCs 1-10
library(data.table)
library(GENESIS)
library(SNPRelate)
library(stringi)
library(readxl)

act <- as.data.table(
  read_xlsx("raw_data/ACT_235_2023_08_21/demo_2023_08_21.xlsx",sheet = 1)
)
act_np <- as.data.table(
  read_xlsx("raw_data/ACT_235_2023_08_21/np_2023_08_21.xlsx", sheet = 1)
)
act_cohorts <- fread(
  "data/act/act_ids_adgc.txt", 
  colClasses = c("character", "character", "factor")
)

load('data/act/mypcair.Rdata')
n_pcs <- 10
pcs <- as.data.table(mypcair$vectors[, 1:(min(32, n_pcs))], 
                     keep.rownames = TRUE)
colnames(pcs) <- stri_replace_all_fixed(colnames(pcs), 'V', 'pc')
setnames(pcs, 'rn', 'IID')

# from https://www.maelstrom-research.org/variable/act_recruitment:gender:Collected
# gender coded as 2 female, 1 male in ACT

act = act[, .(IDFromDave, gender, age_at_death)]
act[gender == 2, gender := 0]
setnames(act, colnames(act), c("IID", "msex", "age_death"))
# act_np[, age_death2 := age_death^2]

act_cohort_dummy <- model.matrix(~ V3, act_cohorts)
act_cohorts <- cbind(act_cohorts, act_cohort_dummy)
act_cohorts <- act_cohorts[, .(V1, V2, V32, V33)]
setnames(act_cohorts, colnames(act_cohorts), 
         c("FID", "IID", "ACT2", "ACT3"))

act <- merge(act, 
                act_cohorts,
                by = "IID")

act <- merge(act, 
                pcs, 
                by = "IID")

act <- act[complete.cases(act)]
setcolorder(act, "FID")
act_np[, FID := 2]

act <- act[!(duplicated(IID))]

act[, FID := "2"] # designated FID for ACT participants

fwrite(act, 
       file = "data/act/act_np.covar",
       sep = " ", 
       quote = FALSE)
