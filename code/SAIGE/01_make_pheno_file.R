#===================================================
# merge .pheno and .covar files for SAIGE analysis
#===================================================

library(data.table)

pheno <- fread("data/mega/mega_np.pheno", na.strings = "-1")
covar <- fread("data/mega/mega_np.covar")

pheno <- merge(pheno, covar)

fwrite(pheno, 
       file = "data/mega/mega_np_pheno_covar.txt", 
       sep = " ",
       na = NA,
       quote = FALSE)
