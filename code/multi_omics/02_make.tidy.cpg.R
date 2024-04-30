##-------------------------------------------
## Pull DNA methylation data from Synapse 
## and prepare for analysis
##-------------------------------------------

library(data.table)
library(GEOquery)
library(synapser)
library(magrittr)
library(readxl)

# login to synapse
synapser::synLogin()

##-----------------------------------------------------
## Connect the specimenIDs in methylation files
## to the projids in our ROSMAP clinical and NP files
##-----------------------------------------------------

# get DNA methylation covariates for each participant (includes batch)
# specimenID
covars <- fread(synapser::synGet("syn5843544")$path, 
                integer64 = "double")
setnames(covars, "Sample", "specimenID")

# get metadata for all biospeciments
# specimenID and individualID
biospecimen.metadata <- fread(synGet("syn21323366")$path)

# create subset for only methylation samples
cpg.iid.metadata <- biospecimen.metadata[notes == "chromatinActivity (methylationArray)"]

# get ROSMAP clinical file from synapse
# individualID & projid
clinical <- fread(synGet("syn3191087")$path)

# merge covars and cpg.iid.metadata to connect specimenID -> individualID
covars.cpg.iid.metadata.tmp <- merge(covars, cpg.iid.metadata, "specimenID")

# merge that data with clinical to connect individualID -> projid
cpg.iid.data <- merge(covars.cpg.iid.metadata.tmp, clinical[, .(projid, individualID)], "individualID")
setnames(cpg.iid.data, "projid", "IID")
cpg.iid.data[, IID := as.character(IID)]

rm(covars, clinical, biospecimen.metadata, cpg.iid.metadata, covars.cpg.iid.metadata.tmp)

##-----------------------------------------------
## See how many methylation samples we have 
## with 1) B-ASC data and 2) in GWAS dataset
##-----------------------------------------------

# load ROSMAP GWAS pheno dataset
load("/home/commons/arter/ROSMAP/01_data/rosmap.ordinal.RData")

# merge ROSMAP GWAS pheno dataset with cpg set to see how many participants in 
# methylation/GWAS overlap (407)
rosmap.cpg <- merge(rosmap, cpg.iid.data[, .(IID, batch, specimenID)], "IID")
rosmap.cpg[, .N]

# load ROSMAP clinical file on statgen to see how many observations total
# that have B-ASC grading (should be 740 I'd think?)
rosmap1 <- as.data.table(read_xlsx("/data_global/ROSMAP/greg_20200109/dataset_843_basic_01-09-2020.xlsx"))
setnames(rosmap1, "projid", "IID")

# merge with cpg.iid.data
# N = 706 with B-ASC grading, 641 >= 80 at death
rosmap1.cpg <- merge(rosmap1, cpg.iid.data[, .(IID, batch, specimenID)], "IID")
rosmap1.cpg[!is.na(arteriol_scler) & age_death > 79, .N]

rm(rosmap, rosmap1)

##---------------------------------------------
## get methylation data from synapse and tidy
## for analysis
##---------------------------------------------

# load methylation data for each cpg site (large .gz file)
if (!file.exists("01_data/ROSMAP_arrayMethylation_imputed.tsv")) {
  gunzip(synGet("syn3168763")$path, 
         destname = "/home/commons/arter/DNA.methylation/01_data/ROSMAP_arrayMethylation_imputed.tsv")
}

cpg <- fread("01_data/ROSMAP_arrayMethylation_imputed.tsv") %>% 
  t()

cpg.colnames <- cpg[1, ]
cpg.specimenIDs <- rownames(cpg)
cpg <- cpg[-1, ] %>% 
  apply(2, as.numeric) %>% 
  as.data.table()
cpg[, specimenID := cpg.specimenIDs[-1]]
setcolorder(cpg, "specimenID")
colnames(cpg)[2:ncol(cpg)] <- cpg.colnames

rm(cpg.colnames, cpg.specimenIDs)

##--------------------------------------------------------------
## Add projids (renamed as IID), batch number, and
## bisulfite conversion efficiency to tidy methylation dataset
##--------------------------------------------------------------

# projids and batch number
cpg <- merge(cpg, cpg.iid.data[, .(IID, specimenID, batch)], "specimenID")

# bisulfite conversion efficiency
load("01_data/bisulfite.conversion.efficiency.RData")
cpg <- merge(cpg, bisulfite.conversion.efficiency, by="specimenID")
setcolorder(cpg, c("IID", "specimenID", "batch", "conversion.efficiency"))

# save dataset
save(cpg, file="01_data/cpg.RData")

