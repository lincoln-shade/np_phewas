library(synapser)
library(synapserutils)

synLogin()
# clin_meta <- syncFromSynapse('syn23277389',
#                          path="data/mcbb/",
#                          ifcollision="keep.local")
# 
# plink_bed <- syncFromSynapse('syn8719897',
#                              path="data/mcbb/",
#                              ifcollision="keep.local")
# 
# plink_bim <- syncFromSynapse('syn8719929',
#                              path="data/mcbb/",
#                              ifcollision="keep.local")
# 
# plink_fam <- syncFromSynapse('syn8719930',
#                              path="data/mcbb/",
#                              ifcollision="keep.local")

syncFromSynapse('syn4650265',
                path="data/mcbb/",
                ifcollision="keep.local")