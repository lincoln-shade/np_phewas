#=================================
# fix NACC IIDs
#=================================

pacman::p_load(data.table, magrittr, stringi)

fam <- fread("data/tmp/act_no_dup.fam")

fam[, V2 := stri_replace_first_regex(V2, ".*_", "")]
# file.rename("data/plink/adc.fam", "data/plink/adc_original_iid.fam")
fwrite(fam, "data/tmp/act_no_dup.fam", col.names = FALSE, sep = " ", quote = FALSE)
