#===========================================
# create list of duplicate IIDs to remove
#===========================================

pacman::p_load(data.table, magrittr, stringi)

psam <- fread("tmp/adc_no_dup.psam", check.names = T)

psam[, IID := stri_replace_first_regex(X.IID, "[[:digit:]]*_", "")]
psam[, IID := stri_replace_first_regex(IID, "_.*", "")]
dup_ids <- psam[duplicated(IID), .(X.IID, SEX)]

if (nrow(dup_ids) > 0) {
  fwrite(dup_ids,
         file = "tmp/adc_dup_naccid.tmp",
         quote = FALSE,
         col.names = FALSE,
         sep = " ")
}
