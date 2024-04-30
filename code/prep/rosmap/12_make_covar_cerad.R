library(data.table)

covar = fread("data/rosmap/rosmap_np.covar")
pheno = fread("data/rosmap/rosmap_np.pheno", na.strings = "-1")

covar = merge(covar, pheno[, .(IID, cerad)], "IID")
covar = covar[!is.na(cerad)]
fwrite(covar, file = "data/rosmap/rosmap_np_cerad.covar")
