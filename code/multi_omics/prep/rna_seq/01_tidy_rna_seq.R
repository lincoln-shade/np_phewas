##-------------------------------------------
## Pull mRNA expression data from Synapse 
## and prepare for analysis
##-------------------------------------------

library(data.table)
library(magrittr)
library(stringi)

##-----------------------------------------------------
## Connect the specimenIDs in RNA expression files
## to the projids in our ROSMAP clinical and NP files
##-----------------------------------------------------

# RNA expression covariate file
covar <- fread("data/synapse/rna_seq/ROSMAP_assay_rnaSeq_metadata.csv")

# get metadata for all biospeciments
# specimenID and individualID
biospec_md <- fread("data/synapse/rna_seq/ROSMAP_biospecimen_metadata.csv")

# create subset of biospecimen metadata for bulk RNAseq in DLPFC
rna_seq_md <- biospec_md[(assay == "rnaSeq") & (nucleicAcidSource == "bulk cell")]

# get ROSMAP clinical file from synapse
# individualID & projid
clinical <- fread("data/synapse/ROSMAP_clinical.csv")

# merge covar with RNAseq metadata
# specimenID < - > individualID
covar_rna_seq_md <- merge(covar, rna_seq_md, "specimenID")

# merge resulting file with the clinical file
# individualID < - > projid
rna_seq_id <- merge(covar_rna_seq_md, clinical[, .(individualID, projid)], "individualID")
rm(covar, clinical, biospec_md, rna_seq_md, covar_rna_seq_md)
# setnames(rna_seq_id, "projid", "IID")

# QC
# keep RNA integrity number > 5
rna_seq_id <- rna_seq_id[RIN > 5]
# what to do about duplicated IIDs??
# seems like main difference is libraryPrep
# polyAselectin vs. rRNAdepletion

##---------------------------------------------
## get normalized mRNA expression
## data from synapse and tidy for analysis
##---------------------------------------------
# rnaseq <- fread(synGet("syn4009614")$path)
rna_seq_1_6 <- 
    fread("data/synapse/rna_seq/ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv") %>% 
    t()

rna_seq_7_8 <- 
    fread("data/synapse/rna_seq/ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv") %>% 
    t()

# test to see if tracking_id is the same as gene_id and
# if both plates 1-6 and 7-8 share the same columns.
# They are, can rbind into one matrix and just keep one of the gene_id rows
ncol(rna_seq_1_6) - sum(as.vector(rna_seq_1_6[1, ]) == as.vector(rna_seq_1_6[2, ]))
ncol(rna_seq_7_8) - sum(as.vector(rna_seq_7_8[2, ]) == as.vector(rna_seq_1_6[2, ]))
rna_seq <- rbind(rna_seq_1_6[-1, ], rna_seq_7_8[-c(1, 2), ])
rm(rna_seq_1_6, rna_seq_7_8)

make_tidy <- function(m) {
    gene_id <- as.vector(m[1, ])
    specimenID <- rownames(m)[-1]
    m <- m[-1, ] %>% 
        apply(2, as.numeric) %>% 
        as.data.table()
    m[, specimenID := specimenID]
    setcolorder(m, "specimenID")
    
    colnames(m)[2:ncol(m)] <- gene_id
    m
}

rna_seq <- make_tidy(rna_seq)

# Can't merge without removing appended batch numbers at the end of the 
# specimenID in the rna_seq data set
rna_seq[, specimenID := stri_replace_last_regex(specimenID, "_[:digit:]*", "")]
rna_seq <- merge(rna_seq,
                 rna_seq_id[, .(specimenID, projid, RIN, libraryPrep)], 
                 "specimenID")
setcolorder(rna_seq, c("projid", "specimenID", "RIN", "libraryPrep"))

# save rna.seq
out_dir <- "data/rna_seq/"
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}
saveRDS(rna_seq, file = "data/rna_seq/rna_seq.Rds")
fwrite(rna_seq, file = "data/rna_seq/rna_seq.csv")
