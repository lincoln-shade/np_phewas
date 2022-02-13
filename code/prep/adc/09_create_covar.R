#==============================================================
# Create covar file
#==============================================================
pacman::p_load(data.table, GENESIS, SNPRelate, stringi)
load('data/adc/mypcair.Rdata')
n_pcs <- 10
pcs <- as.data.table(mypcair$vectors[, 1:(min(32, n_pcs))], 
                     keep.rownames = TRUE)
colnames(pcs) <- stri_replace_all_fixed(colnames(pcs), 'V', 'pc')
setnames(pcs, 'rn', 'IID')

np <- readRDS("data/adc/np.Rds")
# variables to keep
vars <- c("FID", "IID", "NPSEX", "NACCDAGE")

np <- np[, ..vars]
np[, msex := ifelse(NPSEX == 2, 0, NPSEX)]
np[, age_death2 := NACCDAGE^2]
setnames(np, 'NACCDAGE', 'age_death')
np <- merge(np[, .(FID, IID, msex, age_death, age_death2)], pcs, by = "IID")
setcolorder(np, "FID")
##--------------------------------------------
## add indicator variables for ADGC cohorts 
##--------------------------------------------

# load ADGC cohort file
adc <- fread("data/adc/adc_ids_adgc.txt")

np <- merge(np, adc[, .(V2, V3)], by.x = "IID", by.y = "V2") # merge data
setnames(np, "V3", "adgc_adc")
# remove duplicate IIDs, keeping whichever row is in later cohort 
setorder(np, -adgc_adc, IID)
np <- np[!duplicated(IID)]

adgc_adc <- table(np$adgc_adc, useNA = "a")
adgc_adc #check to see no NA
adgc_adc <- adgc_adc[-length(adgc_adc)] # remove NA index
# remove the largest adgc_adc group so that is is the reference group
adgc_adc <- adgc_adc[-(which(adgc_adc == max(adgc_adc)))] 

adgc_adc.names <- names(adgc_adc)

for (i in 1:length(adgc_adc.names)) {
  cohort.num <- paste0("adgc_", adgc_adc.names[i])
  np[, paste(cohort.num) := ifelse(adgc_adc == adgc_adc.names[i], 1, 0)]
  if (adgc_adc[i] != table(np[, ..cohort.num])[2]) {print(paste0(
    "error: cohort ", adgc_adc.names[i], ": ", adgc_adc[i], "; ", cohort.num, ": ", table(np[, ..cohort.num])[2]
  ))}
}

np <- np[, -"adgc_adc"]
setcolorder(np, "FID")

fwrite(np, file = "data/adc/adc_np.covar", sep = " ", quote = FALSE)

