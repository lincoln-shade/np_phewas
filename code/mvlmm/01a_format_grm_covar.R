library(data.table)

grm <- fread("tmp/mega_mvlmm_adnc.grm.gz")
grm_ids <- fread("tmp/mega_mvlmm_adnc.grm.id")
grm_ids[, I := .I]
setnames(grm_ids, "V2", "ID")

grm <- merge(grm, grm_ids[, I, ID], by.x = "V1", by.y = "I")
setnames(grm, "ID", "ID1")

grm <- merge(grm, grm_ids[, I, ID], by.x = "V2", by.y = "I")
setnames(grm, "ID", "ID2")

grm <- grm[, .(ID1, ID2, V4)]

fwrite(grm, file = "tmp/mega_mvlmm_adnc_grm.txt", 
       quote = FALSE, col.names = FALSE, sep = ' ')

pheno <- fread("tmp/mega_mvlmm_adnc.fam")
covar <- fread("data/mega/mega_np.covar")
covar <- merge(pheno[, .(V2)], covar, by.x = "V2", by.y = "IID")
covar[, int := 1]
setcolorder(covar, 'int')

fwrite(covar, file = "tmp/mega_mvlmm_adnc.covar", 
       sep = ' ', quote = FALSE, col.names = FALSE)
