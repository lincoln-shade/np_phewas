#===========================================
# create list of duplicate IIDs to remove
#===========================================

pacman::p_load(data.table, magrittr, stringi)

fam <- fread("tmp/adc_no_naccid_dups.fam")

naccids <- fam[, V2 := stri_replace_first_regex(V2, ".*_", "")]

fwrite(naccids, 
       file = "tmp/adc_no_naccid_dups.fam",
       col.names = FALSE,
       sep = " ")
