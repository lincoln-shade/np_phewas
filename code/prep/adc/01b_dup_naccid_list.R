#===========================================
# create list of duplicate IIDs to remove
#===========================================

pacman::p_load(data.table, magrittr, stringi)

fam <- fread("tmp/adc_no_dup.fam")

fam[, IID := stri_replace_first_regex(V2, ".*_", "")]
dup_ids <- fam[duplicated(IID), .(V1, V2)]

if (nrow(dup_ids) > 0) {
  fwrite(dup_ids,
         file = "tmp/adc_dup_naccid.tmp",
         quote = FALSE,
         col.names = FALSE,
         sep = " ")
}