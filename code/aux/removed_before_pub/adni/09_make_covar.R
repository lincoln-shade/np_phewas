#====================================================
# create covariate file
#====================================================

library(data.table)
library(magrittr)
library(stringi)

np <- fread("raw_data/NEUROPATH_07_06_21.csv")
fam <- fread("data/adni/adni_np.fam")
demo <- fread("raw_data/PTDEMOG.csv") 
demo <- demo[!duplicated(RID)]
np <- merge(fam, np, by.x = 'V2', by.y = 'RID')
np <- merge(np, demo, by.x = 'V2', by.y = 'RID')
load('data/adni/mypcair.Rdata')

n_pcs <- 10
pcs <- as.data.table(mypcair$vectors[, 1:(min(32, n_pcs))], 
                     keep.rownames = TRUE)
colnames(pcs) <- stri_replace_all_fixed(colnames(pcs), 'V', 'pc')
setnames(pcs, 'rn', 'IID')

# variables to keep
setnames(np, c('V1', 'V2', 'PTGENDER', 'NPDAGE'), c('FID', 'IID', 'msex', 'age_death'))
vars <- c("FID", "IID", "msex", "age_death")
np[, IID := as.character(IID)]
np <- np[, ..vars]
np[, msex := ifelse(msex == 2, 0, msex)]
np[, age_death2 := age_death^2]
np <- merge(np, pcs, by = "IID")
setcolorder(np, "FID")

fwrite(np, file = "data/adni/adni_np.covar", sep = " ", quote = FALSE)
