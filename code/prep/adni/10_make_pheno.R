#====================================================
# create covariate file
#====================================================

library(data.table)
library(magrittr)
library(stringi)

np <- fread("raw_data/NEUROPATH_07_06_21.csv")
fam <- fread("data/adni/adni_np.fam")
np <- merge(fam, np, by.x = 'V2', by.y = 'RID')

setnames(np, c('V1', 'V2'), c('FID', 'IID'))

# turn lewy body variable from categorical to binary
# NPLBOD 0 | 1 2 3 4
# NPAVAS 0 | 1 2 3, 9 --> NA, only 6 controls
# NPOLD microinfarcts 0 | 1
# NPINF gross infarcts 0 | 1, only 4 cases
# NPHISCL 0 | 1 3, only 6 cases
# NPBRAAK 0 1 2 3 4 | 5 6, 13 controls
# NPNEUR 0 1 2 | 3, 16 controls
# NPAMY CAA 0 | 1 2 3, only 8 controls
# NPARTER B-ASC 0 1 | 2 3, 14 cases

np[NPLBOD == 0, lewy := 0]
np[NPLBOD %in% 1:4, lewy := 1]

np[NPAVAS == 0, athero := 0]
np[NPAVAS %in% 1:3, athero := 1]

np[NPINF == 0, grossinf := 0]
np[NPINF %in% 1, grossinf := 1]

np[NPOLD == 0, microinf := 0]
np[NPOLD %in% 1, microinf := 1]

np[NPHIPSCL == 0, hs := 0]
np[NPHIPSCL %in% 1:3, hs := 1]

np[NPBRAAK %in% 0:4, braak56 := 0]
np[NPBRAAK %in% 5:6, braak56 := 1]

np[NPNEUR %in% 0:2, cerad3 := 0]
np[NPNEUR %in% 3, cerad3 := 1]

np[NPAMY %in% 0, caa := 0]
np[NPAMY %in% 1:3, caa := 1]

np[NPARTER %in% 0:1, arteriol23 := 0]
np[NPARTER %in% 2:36, arteriol23 := 1]

np_vars <- c('FID', 'IID', 'lewy', 'athero', 'grossinf', 'microinf', 'hs', 
             'braak56', 'cerad3', 'caa', 'arteriol23')
np <- np[, ..np_vars]

fwrite(np, 
       file = "data/adni/adni_np.pheno", 
       sep = " ", 
       quote = FALSE, 
       na = '-1')
