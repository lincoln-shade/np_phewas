# merge phenotype files for NACC, ROSMAP, and ACT
library(data.table)

nacc_ord = fread("shared/adc_np_ord.txt", na.strings = "-1")
nacc_bin = fread("shared/adc_np_bin.txt", na.strings = "-1")
rosmap = fread("shared/rosmap_np.pheno", na.strings = "-1")
act = fread("data/act/act_np.pheno", na.strings = "-1")

# merge nacc_ord and nacc_bin
nacc = merge(nacc_ord, nacc_bin, c("FID", "IID"), all.x = TRUE, all.y = TRUE)

# merge phenotype files
pheno = rbindlist(list(nacc, rosmap, act), use.names = TRUE)

fwrite(pheno, file = "shared/nacc_rosmap_act_np.pheno")
