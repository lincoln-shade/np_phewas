library(data.table)

covar = fread("data/act/act_np.covar")
pheno = fread("data/act/act_np.pheno", na.strings = "-1")

covar = merge(covar, pheno[, .(IID, cerad)], "IID")
covar = covar[!is.na(cerad)]
fwrite(covar, file = "data/act/act_np_cerad.covar")
