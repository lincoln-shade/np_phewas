library(data.table)
library(synapser)
library(magrittr)
library(stringi)

# array metadata: syn23569440

# login to synapse
synapser::synLogin()

# RNA expression metadata file
rnaseq_meta <- fread(synGet("syn21088596")$path)
rnaseq_meta <- rnaseq_meta[!is.na(RIN)][libraryPrep == "polyAselection"]

# get metadata for all biospeciments
# specimenID and individualID
biospec_meta <- fread(synGet("syn21323366")$path)
biospec_meta <- biospec_meta[(assay == "rnaSeq") & 
                               (nucleicAcidSource == "bulk cell") &
                               (tissue == "dorsolateral prefrontal cortex")]

# get ROSMAP clinical file from synapse
# individualID & projid
clinical_meta <- fread(synGet("syn3191087")$path)

# merge covars with RNAseq metadata
# specimenID < - > individualID
meta <- merge(rnaseq_meta, biospec_meta, "specimenID")
meta <- merge(meta, clinical_meta, "individualID")

##---------------------------------------------
## get normalized mRNA expression
## data from synapse and tidy for analysis
##---------------------------------------------

rnaseq_exp1_6 <- fread(synGet("syn3505732")$path)
rnaseq_exp7_8 <- fread(synGet("syn3505724")$path) 

# remove label columns and transpose
tidy_rna_exp <- function(dt, meta_vars=c("tracking_id", "gene_id")) {
  dt_meta <- dt[, ..meta_vars]
  dt_data <- dt[, -c(..meta_vars)] %>% 
    t() %>% 
    as.data.table(keep.rownames = "IID")
  colnames(dt_data)[2:ncol(dt_data)] <- dt_meta[[1]]
  return(dt_data)
}

rnaseq_exp1_6_tidy <- tidy_rna_exp(rnaseq_exp7_8)
rnaseq_exp7_8_tidy <- tidy_rna_exp(rnaseq_exp1_6)

rnaseq_tidy <- rbind(rnaseq_exp1_6_tidy, rnaseq_exp7_8_tidy)
rm(rnaseq_exp1_6, rnaseq_exp1_6_tidy, rnaseq_exp7_8, rnaseq_exp7_8_tidy)

# convert specimenID to projid and add RIN and PMI covariates
rnaseq_tidy[, IID := stri_replace_last_regex(IID, "_[0-9]$", "")]
meta <- meta[specimenID %in% rnaseq_tidy$IID]
rnaseq_tidy <- rnaseq_tidy[IID %in% meta$specimenID]
rnaseq_tidy <- rnaseq_tidy[!duplicated(IID)]
rnaseq_tidy <- merge(meta[, .(specimenID, projid, RIN, pmi, msex, age_death, 
                              braaksc)], 
                     rnaseq_tidy, 
                     by.x = "specimenID",
                     by.y = "IID")
rnaseq_tidy[, specimenID := NULL]

fwrite(rnaseq_tidy, file = "data/rnaseq/rnaseq.txt.gz", 
       quote = FALSE, 
       sep = " ",
       na = "NA")
