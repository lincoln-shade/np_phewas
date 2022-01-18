#==============================================================
# perform PC-AiR and create covar file
#==============================================================


pacman::p_load(data.table, GENESIS, SNPRelate)

#-----------
# PC-AiR
#-----------
# mostly taken from here: 
# https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html
gdsfile <- "tmp/pcair.gds"
plink_path <- "data/adc/adc_np_pruned"
snps <- fread(paste0(plink_path, ".bim"))
snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"),
              bim.fn = paste0(plink_path, ".bim"),
              fam.fn = paste0(plink_path, ".fam"),
              out.gdsfn = gdsfile)

# create kinship matrix
KINGmat <- kingToMatrix("data/adc/adc_np.kin", estimator = "Kinship")

geno <- GWASTools::GdsGenotypeReader(filename = gdsfile)

mypcair <- pcair(geno, kinobj = KINGmat, divobj = KINGmat, snp.include = snps$V2)

pcs <- as.data.table(mypcair$vectors[, 1:5], keep.rownames = TRUE)
setnames(pcs, colnames(pcs), c("IID", "PC1", "PC2", "PC3", "PC4", "PC5"))

#---------------------------------
# create pheno analysis data set
#---------------------------------

np <- readRDS("data/adc/np.Rds")
# variables to keep
vars <- c("FID", "IID", "NPSEX", "NACCDAGE")

np <- np[, ..vars]
np <- merge(np, pcs[, 1:6], by = "IID")
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
##--------------------
## write pheno files
##--------------------
# ordinal
fwrite(np, file = "data/adc/adc_np.covar", sep = " ", quote = FALSE)

