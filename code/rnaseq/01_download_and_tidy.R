##-------------------------------------------
## Pull mRNA expression data from Synapse 
## and prepare for analysis
##-------------------------------------------

library(data.table)
library(synapser)
library(magrittr)
library(stringi)

# array metadata: syn23569440

# login to synapse
synapser::synLogin()

# download mRNA array data from Synapse
rna_array <- fread(synGet("syn4009614")$path)

# RNAseq metadata file with specimenID (RIN missing in RNA array metadata file)
rnaseq_meta <- fread(synGet("syn21088596")$path)
rnaseq_meta <- rnaseq_meta[!is.na(RIN)]
# biospecimen metadata file with individualID and specimenID
biospecimen_meta <- fread(synGet("syn21323366")$path) 
biospecimen_meta <- biospecimen_meta[assay == "rnaSeq", 
                                     .(individualID, specimenID)]
biospecimen_meta <- biospecimen_meta[!duplicated(biospecimen_meta)]
meta <- merge(rnaseq_meta, biospecimen_meta, by = "specimenID")
# get clincal data with projid and individualID
clinical <- fread(synGet("syn3191087")$path)
clinical[, IID := as.character(projid)]

meta <- merge(meta, clinical, by = "individualID")

# tidy rna_array
ids <- colnames(rna_array)
rna_array_meta <- rna_array[, 1:3]
rna_array_data <- rna_array[, 4:ncol(rna_array)]
rna_arrayt <- t(rna_array_data)
rna_arrayt_dt <- as.data.table(rna_arrayt, keep.rownames = "IID")
colnames(rna_arrayt_dt)[2:ncol(rna_arrayt_dt)] <- 
  rna_array_meta[, paste0(Symbol, "_", ProbeID)]

rna_arrayt_dt[, IID := stri_replace_first_regex(IID, ".*_", "")]

