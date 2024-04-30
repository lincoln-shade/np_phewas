
# dichotomize ordinal phenotypes for colocalization analysis
library(data.table)

phenos <- fread("shared/nacc_rosmap_act_np.pheno", na.strings = "")

# braak stage: 0-3 vs 4-6
phenos[braak < 4, braak := 0]
phenos[braak > 0, braak := 1]

# cerad score: 0-1 vs 2-3
phenos[cerad < 2, cerad := 0]
phenos[cerad > 0, cerad := 1]

# # diffuse_abeta: 0-1 vs 2-3
phenos[diffuse_abeta < 2, diffuse_abeta := 0]
phenos[diffuse_abeta > 0, diffuse_abeta := 1]

# arteriol 0-1 vs 2-3
phenos[arteriol < 2, arteriol := 0]
phenos[arteriol > 0, arteriol := 1]

# atheroscler: 0-1 vs 2-3
phenos[atheroscler < 2, atheroscler := 0]
phenos[atheroscler > 0, atheroscler := 1]

# late: 0-1 vs 2-3
phenos[late < 2, late := 0]
phenos[late > 0, late := 1]

# lewy_body: 0-2 vs 3
phenos[lewy_body < 3, lewy_body := 0]
phenos[lewy_body > 0, lewy_body := 1]

# caa_ord: 0-1 vs 2-3
phenos[caa < 2, caa := 0]
phenos[caa > 0, caa := 1]

fwrite(phenos, "shared/nacc_rosmap_act_np_dichotomized.pheno")
