#===========================================
# create list of duplicate IIDs to remove
#===========================================

pacman::p_load(data.table, magrittr, stringi)

psam <- fread("tmp/adc_no_naccid_dups.psam")

naccids <- psam[, `#IID` := 
                  stri_replace_first_regex(`#IID`, ".*_", "")]
naccids <- psam[, `#IID` := stri_replace_first_regex(`#IID`, "_.*", "")]
naccids[, `#FID` := "0"]
setnames(naccids, "#IID", "IID")
setcolorder(naccids, c("#FID", "IID", "SEX"))
fwrite(naccids, 
       file = "tmp/adc_no_naccid_dups.psam",
       sep = " ",
       na = 'NA',
       quote = FALSE)
