#==============================================================
# perform PC-AiR and create covar file
#==============================================================

library(pacman)
p_load(data.table, GENESIS, SNPRelate)

#-----------
# PC-AiR
#-----------
# (mostly taken from here: https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html)
gdsfile <- "data/tmp/pcair.gds"
plink_path <- "data/plink/adc_np_pruned"
snps <- fread(paste0(plink_path, ".bim"))
snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"),
              bim.fn = paste0(plink_path, ".bim"),
              fam.fn = paste0(plink_path, ".fam"),
              out.gdsfn = gdsfile)

# # create list of uncorrelated SNPs
# 
# gds <- snpgdsOpen(gdsfile)
# snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
#                           ld.threshold=sqrt(0.1), verbose = TRUE)
# pruned <- unlist(snpset, use.names=FALSE)
# snpgdsClose(gds)

# create kinship matric
KINGmat <- kingToMatrix("data/adc_np.kin", estimator = "Kinship")

geno <- GWASTools::GdsGenotypeReader(filename = gdsfile)

mypcair <- pcair(geno, kinobj = KINGmat, divobj = KINGmat, snp.include = snps$V2)

pcs <- as.data.table(mypcair$vectors[, 1:5], keep.rownames = TRUE)
setnames(pcs, colnames(pcs), c("IID", "PC1", "PC2", "PC3", "PC4", "PC5"))

#---------------------------------
# create pheno analysis data set
#---------------------------------

np <- readRDS("data/np.Rds")
# variables to keep
vars <- c("FID", "IID", "NPSEX", "npDAGE")

np <- np[, ..vars]
np <- merge(np, pcs[, 1:6], by = "IID")
setcolorder(np, "FID")
##--------------------------------------------
## add indicator variables for ADGC cohorts 
##--------------------------------------------

# load ADGC cohort file
adc <- fread("data/adc_ids_adgc.txt")

np <- merge(np, adc[, .(V2, V3)], by.x = "IID", by.y = "V2") # merge data
setnames(np, "V3", "adgc_adc")
# remove duplicate IIDs, keeping whichever row is in later cohort 
setorder(np, -adgc_adc, IID)
np <- np[!duplicated(IID)]

adgc_adc <- table(np$adgc_adc, useNA = "a")
adgc_adc #check to see no NA
adgc_adc <- adgc_adc[-length(adgc_adc)] # remove NA index
adgc_adc <- adgc_adc[-(which(adgc_adc == max(adgc_adc)))] # remove the largest adgc_adc group so that is is the "control" group without an indicator

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
##--------------------
## write pheno files
##--------------------
# ordinal
fwrite(np, file = "data/plink/adc_np.covar", sep = " ", quote = FALSE)

rm(list = ls())
p_unload(all)
