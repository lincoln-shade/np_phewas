library(data.table)

covar = fread("data/adc/adc_np.covar")
pheno = fread("data/adc/adc_np_ord.txt", na.strings = "-1")

covar = merge(covar, pheno[, .(IID, cerad)], "IID")
covar = covar[!is.na(cerad)]
covar[, cerad1 := ifelse(cerad == 1, 1, 0)]
covar[, cerad2 := ifelse(cerad == 2, 1, 0)]
covar[, cerad3 := ifelse(cerad == 3, 1, 0)]
covar[, cerad := NULL]
fwrite(covar, file = "data/adc/adc_np_cerad.covar")
