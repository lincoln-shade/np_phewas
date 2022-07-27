#===========================================
# create list of duplicate IIDs to remove
#===========================================

pacman::p_load(data.table, magrittr, stringi)

psam <- fread("tmp/adc_no_naccid_dups.psam")

naccids <- psam[, `#IID` := 
                  stri_replace_first_regex(`#IID`, "[[:digit:]]*_", "")]
naccids <- psam[, `#IID` := stri_replace_first_regex(`#IID`, "_.*", "")]
fwrite(naccids, 
       file = "tmp/adc_no_naccid_dups.psam",
       sep = " ",
       na = 'NA',
       quote = FALSE)
