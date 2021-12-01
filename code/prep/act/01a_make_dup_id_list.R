#===============================================
# Create list of duplicates to remove
#===============================================

library(pacman)
p_load(data.table, magrittr, stringi)
cargs <- commandArgs(trailingOnly = TRUE)
prefix <- "act"
np_ids <- fread("data/act_np_ids.txt", header = FALSE)
miss <- fread(paste0("data/tmp/", prefix, "_dup.imiss"))
fam <- fread("data/plink/act.fam", header = FALSE)
related <- fread(paste0("data/tmp/", prefix, "_dup.con"))
setnames(related, c("ID1", "ID2"), c("IID1", "IID2"))
related[, PI_HAT := (2 * N_IBS2 + N_IBS1) / (2 * (N_IBS0 + N_IBS1 + N_IBS2))]
related[, pair := 1:.N]
related <- related[, .(IID1, IID2, pair, PI_HAT)]
related[, IID1_short := stri_replace_first_regex(IID1, ".*_", "")]
related[, IID2_short := stri_replace_first_regex(IID2, ".*_", "")]

related_long <- melt.data.table(related, 
                                measure.vars = c("IID1", "IID2"), 
                                value.name = "IID")
related_long <- merge(related_long, miss[, .(IID, F_MISS)], "IID")
# identify ids in multiple pairs
dup_iids <- related_long$IID[duplicated(related_long$IID)]
setorder(related_long, pair, variable)

# related pairs with uniq "short" ids 
# (i.e. one id not the other after .*_ prefix removed)
uniq_id_pairs <- related[IID1_short != IID2_short, pair]
related_long_uniq <- related_long[pair %in% uniq_id_pairs]
related_long_dup <- related_long[!(pair %in% uniq_id_pairs)]

#----------------------------------
# pairs with uniq iids
# prefer iids in neuropath data set
# (only 1 pair of this type)
#----------------------------------

np_ids_grep_pattern <- paste(np_ids$V2, collapse = "|")

iid_remove <- related_long_uniq[!(grep(np_ids_grep_pattern, IID)), IID]
#---------------------
# pairs with dup iids
# check missingness
#---------------------
setorder(related_long_dup, pair, F_MISS)
iid_remove <- c(iid_remove, 
                related_long_dup[duplicated(pair), IID])

id_remove <- fam[V2 %in% iid_remove, .(V1, V2)]

write.table(id_remove, 
            file = paste0("data/tmp/", prefix, "_dup_remove.tmp"), 
            quote = F, 
            row.names = F, 
            col.names = F)

rm(list = ls())
p_unload(all)
