library(data.table)
library(magrittr)
source("code/functions/strip_alleles.R")

# rosmap_phenos <- c('arteriol_scler_bin', '')
rosmap_path <- "output/gwas/rosmap"
adc_path <- "output/gwas/adc"

get_file_name <- function(path, pheno) {
  file_name <- paste0(path, '/', pheno, '.assoc.logistic')
}

format_results <- function(file_name) {
  results <- fread(file_name)
  results <- results %>% 
    .[!(is.na(P))] %>% 
    setorder(P)
}

# ASC
asc_ros <- format_results(get_file_name(rosmap_path, 'arteriol_scler_bin'))

asc_adc <- format_results(get_file_name(adc_path, 'ASC'))
asc_adc[, SNP := strip_alleles(SNP)]

asc <- merge(asc_ros, asc_adc, c('CHR', 'BP', 'A1'), suffixes = c('_ros', '_adc'))

asc[P_ros < 0.00005 & P_adc < 0.05, .(CHR, BP, SNP_ros, A1, OR_ros, P_ros, OR_adc, P_adc)]

fwrite(asc_adc, file = "tmp/asc_adc.assoc.logistic", quote = F, sep = ' ', na = 'NA')
asc_meta <- fread("output/asc_meta.meta")
setorder(asc_meta, P)
head(asc_meta)

# microinf
micr_ros <- format_results(get_file_name(rosmap_path, 'ci_num2_mct'))

micr_adc <- format_results(get_file_name(adc_path, 'NACCMICR'))
micr_adc[, SNP := strip_alleles(SNP)]

micr <- merge(micr_ros, micr_adc, c('CHR', 'BP', 'A1'), suffixes = c('_ros', '_adc'))

micr[P_ros < 0.005 & P_adc < 0.0005, .(CHR, BP, SNP_ros, A1, OR_ros, P_ros, OR_adc, P_adc)]

fwrite(micr_adc, file = "tmp/micr_adc.assoc.logistic", quote = F, sep = ' ', na = 'NA')
micr_meta <- fread("output/micr_meta.meta")
setorder(micr_meta, P)
head(micr_meta)

# macroinf
macr_ros <- format_results(get_file_name(rosmap_path, 'ci_num2_gct'))

macr_adc <- format_results(get_file_name(adc_path, 'NACCINF'))
macr_adc[, SNP := strip_alleles(SNP)]

macr <- merge(macr_ros, macr_adc, c('CHR', 'BP', 'A1'), suffixes = c('_ros', '_adc'))

macr[P_ros < 0.005 & P_adc < 0.0005, .(CHR, BP, SNP_ros, A1, OR_ros, P_ros, OR_adc, P_adc)]

fwrite(macr_adc, file = "tmp/macr_adc.assoc.logistic", quote = F, sep = ' ', na = 'NA')
macr_meta <- fread("output/macr_meta.meta")
setorder(macr_meta, P)
head(macr_meta)

# atherosclerosis
ath_ros <- format_results(get_file_name(rosmap_path, 'cvda_4gp2_bin'))

ath_adc <- format_results(get_file_name(adc_path, 'ATHCW'))
ath_adc[, SNP := strip_alleles(SNP)]

ath <- merge(ath_ros, ath_adc, c('CHR', 'BP', 'A1'), suffixes = c('_ros', '_adc'))

ath[P_ros < 0.005 & P_adc < 0.0005, .(CHR, BP, SNP_ros, A1, OR_ros, P_ros, OR_adc, P_adc)]

fwrite(ath_adc, file = "tmp/ath_adc.assoc.logistic", quote = F, sep = ' ', na = 'NA')
ath_meta <- fread("output/ath_meta.meta")
setorder(ath_meta, P)
head(ath_meta)

# hippocampal sclerosis
hs_ros <- format_results(get_file_name(rosmap_path, 'hspath_typ'))

hs_adc <- format_results(get_file_name(adc_path, 'HS'))
hs_adc[, SNP := strip_alleles(SNP)]

hs <- merge(hs_ros, hs_adc, c('CHR', 'BP', 'A1'), suffixes = c('_ros', '_adc'))

hs[CHR != 7][P_ros < 0.005 & P_adc < 0.0005, .(CHR, BP, SNP_ros, A1, OR_ros, P_ros, OR_adc, P_adc)]

fwrite(hs_adc, file = "tmp/hs_adc.assoc.logistic", quote = F, sep = ' ', na = 'NA')
hs_meta <- fread("output/hs_meta.meta")
setorder(hs_meta, P)
head(hs_meta)

# Lewy bodies
lewy_ros <- format_results(get_file_name(rosmap_path, 'dlbdx_bin'))

lewy_adc <- format_results(get_file_name(adc_path, 'LEWY'))
lewy_adc[, SNP := strip_alleles(SNP)]

lewy <- merge(lewy_ros, lewy_adc, c('CHR', 'BP', 'A1'), suffixes = c('_ros', '_adc'))

lewy[P_ros < 0.005 & P_adc < 0.0005, .(CHR, BP, SNP_ros, A1, OR_ros, P_ros, OR_adc, P_adc)]

fwrite(lewy_adc, file = "tmp/lewy_adc.assoc.logistic", quote = F, sep = ' ', na = 'NA')
lewy_meta <- fread("output/lewy_meta.meta")
setorder(lewy_meta, P)
head(lewy_meta)
