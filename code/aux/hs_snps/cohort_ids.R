library(data.table)
library(magrittr)

pheno <- fread('data/mega/mega_np.pheno', na.strings = '-1')
covar <- fread('data/mega/mega_np.covar')
dt <- merge(pheno, covar, c('FID', 'IID'))
dt <- dt[!is.na(hs)]

na_strings <- c("-9", "-4", "9999", "999", "888")
uds <- fread("/data_global/nacc/investigator_nacc56.csv", 
             na.strings = na_strings)

# NP data not longitudinal, so just keep obs from last visit
uds <- uds[NACCVNUM == NACCAVST][, .(NACCID, NPFORMVER, NACCYOD)]
setnames(uds, 'NACCYOD', 'DEATHYR')

# Minimal Data Set + NP
mds <- fread("/data_global/nacc/fardo09062019.csv", 
             na.strings = na_strings)
mds <- mds[, .(NACCID, NPFORMVER, DEATHYR)]

nacc <- rbind(uds, mds)
nacc <- nacc[!is.na(NPFORMVER)]

nacc <- merge(dt, nacc, by.x = 'IID', by.y = 'NACCID')
nacc_old <- nacc[(DEATHYR >= 2000 & DEATHYR < 2014) & 
                   ((adgc_7 == 0) & (adgc_8 == 0) & (adgc_9 == 0) & 
                      (adgc_10 == 0) & (adgc_11 == 0) & (adgc_12 == 0))]
nacc_new <- nacc[!(IID %in% nacc_old$IID) & (DEATHYR >= 2000)]
nacc_b2000_NA <- nacc[DEATHYR < 2000 | is.na(DEATHYR)]

rosmap <- dt[FID == '1' & !is.na(hs)]
act <- dt[FID == '2' & !is.na(hs)]
adni <- dt[FID == 'ADNI' & !is.na(hs)]
nacc_b2000_NA_adni <- rbind(nacc_b2000_NA, adni, fill=TRUE)

fwrite_ids <- function(dt, out_prefix='data/aux/hs_snps/') {
  dt_name <- deparse(substitute(dt))
  fwrite(dt[, .(FID, IID)],
         file = paste0(out_prefix, dt_name, '_ids.txt'),
         quote = FALSE, col.names = FALSE, sep = ' ')
}

fwrite_ids(nacc_old)
fwrite_ids(nacc_new)
fwrite_ids(nacc_b2000_NA_adni)
fwrite_ids(rosmap)
fwrite_ids(act)
