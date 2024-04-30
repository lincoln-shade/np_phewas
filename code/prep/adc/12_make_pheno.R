# create NACC phenotyp file for POLMM analysis
library(data.table)
np = readRDS("data/adc/np_qced.Rds")



# TDP-43 ordinal var
# Stage 1: Amygdala only
# Stage 2: + Hippocampus
# Stage 3: + Middle frontal gyrus or cortex
np[NPTDPB == 0 | NPTDPC == 0 | NPTDPD == 0 | NPTDPE == 0, late := 0]
np[NPTDPB == 1, late := 1] 
np[NPTDPC == 1 | NPTDPD == 1, late := 2]
np[NPTDPE == 1, late := 3] 

# neuritic plaque CERAD score ordinal
np[!is.na(NACCNEUR), cerad := NACCNEUR]

# braak score ordinal
np[!is.na(NACCBRAA), braak := NACCBRAA]

# Thal CERAD score / diffuse plaques ordinal and binary vars
np[!is.na(NACCDIFF), diffuse_abeta := NACCDIFF]

# arteriolosclerosis
np[!is.na(NACCARTE), arteriol := NACCARTE]

# atherosclerosis in circle of willis
np[NACCAVAS %in% 0:3, atheroscler := NACCAVAS]

# lewy body
np[NACCLEWY %in% 0:3, lewy_body := NACCLEWY]
np[NACCLEWY == 4 & NPLBOD == 5, lewy_body := 0] # olfactory

# cerebral amyloid angiopathy
np[NACCAMY %in% 0:3, caa := NACCAMY]

# Hippocampal sclerosis
np[NPHIPSCL == 0 | NPSCL == 0, hs := 0]
np[NPHIPSCL %in% c(1, 2, 3) | NPSCL == 1, hs := 1]

# microinfarcts
np[NACCMICR == 0, microinf := 0]
np[NACCMICR == 1, microinf := 1]

# gross infarcts
np[NACCINF == 0, grossinf := 0]
np[NACCINF == 1, grossinf := 1]

keep_vars_ord = c(
  "FID", "IID", "late", "braak", "cerad", "diffuse_abeta", "arteriol", 
  "atheroscler", "lewy_body", "caa"
)

keep_vars_bin = c("FID", "IID", "hs", "microinf", "grossinf")
np_pheno_ord = np[, ..keep_vars_ord]
np_pheno_bin = np[, ..keep_vars_bin]

# write table
fwrite(np_pheno_ord, file = "data/adc/adc_np_ord.txt", na = "-1")
fwrite(np_pheno_bin, file = "data/adc/adc_np_bin.txt", na = "-1")
