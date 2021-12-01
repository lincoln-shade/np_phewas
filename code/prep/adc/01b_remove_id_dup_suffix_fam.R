#=================================
# fix NACC IIDs
#=================================

pacman::p_load(data.table, magrittr, stringi)

fam <- fread("data/tmp/adc_no_dup.fam")

# one iid (in a threesome) duplicate made it through, so leave as duplicate
fam[V2 != "NACC438637_dup", V2 := stri_replace_first_regex(V2, "_dup", "")]
# file.rename("data/plink/adc.fam", "data/plink/adc_original_iid.fam")
fwrite(fam, "data/tmp/adc_no_dup.fam", col.names = FALSE, sep = " ", quote = FALSE)
