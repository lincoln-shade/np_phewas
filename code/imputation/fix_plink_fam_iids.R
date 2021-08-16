#=================================
# fix NACC IIDs
#=================================

pacman::p_load(data.table, magrittr, stringi)

fam <- fread("data/plink/adc_original_iid.fam")
fam[, V2 := stri_replace_first_regex(V2, "[:alnum:]*_", "")]
fam[, V2 := stri_replace_first_regex(V2, "[:alnum:]*:", "")]
fam[, V2 := stri_replace_first_fixed(V2, "_2", "")]
fam[duplicated(fam), V2 := paste0(V2, "_dup")]
# file.rename("data/plink/adc.fam", "data/plink/adc_original_iid.fam")
fwrite(fam, "data/plink/adc.fam", col.names = FALSE, sep = " ", quote = FALSE)
