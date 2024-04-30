#=================================
# fix NACC IIDs
#=================================

library(data.table)
library(magrittr)
library(stringi)

fam <- fread("tmp/act_no_dup.fam")

fam[, V2 := stri_replace_first_regex(V2, ".*_", "")]
# file.rename("data/plink/adc.fam", "data/plink/adc_original_iid.fam")
fwrite(fam, "tmp/act_no_dup.fam", col.names = FALSE, sep = " ", quote = FALSE)
