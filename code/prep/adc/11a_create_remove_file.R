#===============================================
# Create list of related individuals to remove
#===============================================

library(pacman)
p_load(data.table, magrittr, stringi)
cl_args <- commandArgs(trailingOnly = TRUE)
dot_genome <- cl_args[1] # .genome file from PLINK
np_path <- cl_args[2] # file with participants
phenotype <- cl_args[3] # phenotype

dir_path <- "data/related_rm"
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

related <- fread(dot_genome)
related[, pair := 1:.N]
ids <- rbind(related[, .(FID1, IID1)] %>% 
               setnames(colnames(.), c("FID", "IID")),
             related[, .(FID2, IID2)] %>% 
               setnames(colnames(.), c("FID", "IID")))
related <- related[, .(IID1, IID2, pair, PI_HAT)]
related.long <- melt.data.table(related, measure.vars = c("IID1", "IID2"), value.name = "IID")

np <- readRDS(np_path)
setnames(np, paste(phenotype), "phenotype")
covar <- fread("data/plink/adc_np.covar")
np <- np[IID %in% covar$IID]
related_pheno <- merge(np[, .(FID, IID, phenotype)], related.long, by = "IID")

related_pheno <- related_pheno[!(is.na(phenotype))][pair %in% pair[duplicated(pair)]]

#-------------------------------------------------------
# add random number seed and create list of ids to remove
#-------------------------------------------------------

if (nrow(related_pheno) > 0) {
  set.seed(20210816)
  related_pheno[, n := sample(1:.N, .N)]
  setorder(related_pheno, pair, n)
  
  id_remove <- related_pheno[seq(1, .N, 2), .(FID, IID)]
} else {
  id_remove <- data.table(FID = "", IID = "")
}

fwrite(id_remove, 
       file = paste0(dir_path, "/", phenotype, ".remove"), 
       col.names = FALSE, 
       quote = FALSE,
       sep = " ")

rm(list = ls())
p_unload(all)
