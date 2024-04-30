library(data.table)
library(ggplot2)
library(coloc)

bel2022 <- 
  fread("/data_global/summary_statistics/Bellenguez/GCST90027158_buildGRCh38.tsv.gz")
bel2022 <- bel2022[!which(duplicated(variant_id))]
hs <- fread("output/gwas/mega/saige/hs_saige_results.txt")
hs_pheno <- fread("data/mega/mega_np.pheno", na.strings = "-1")
late <- fread("output/gwas/mega/polmm/late_polmm_results.txt")
late_pheno <- fread("data/mega/mega_np_ord.pheno", na.strings = "-1")
setnames(hs, c("MarkerID", "p.value"), c("variant_id", "p_value"))
setnames(late, c("SNPID", "pval.spa"), c("variant_id", "p_value"))
bel2022[, Phenotype := "AD_Bellenguez_et_al_2022"]
hs[, Phenotype := "HS"]
late[, Phenotype := "LATE_NC"]
coloc_radius <- 2e5
bel2022_tmem <- bel2022[chromosome == 7
][base_pair_location > (12229967 - coloc_radius) & 
    base_pair_location < (12229967 + coloc_radius)]
bel2022_grn <- bel2022[chromosome == 17
][base_pair_location > (44352876 - coloc_radius) & 
    base_pair_location < (44352876 + coloc_radius)]

hs_tmem <- merge(hs[CHR == 7], 
                 bel2022_tmem[, .(variant_id, base_pair_location)], 
                 by = "variant_id")
late_tmem <- merge(late[chr == 7], 
                 bel2022_tmem[, .(variant_id, base_pair_location)], 
                 by = "variant_id")
bel2022_tmem <- bel2022_tmem[variant_id %in% hs_tmem$variant_id]
tmem <- rbindlist(list(
  bel2022_tmem[, .(Phenotype, variant_id, base_pair_location, p_value)],
  hs_tmem[, .(Phenotype, variant_id, base_pair_location, p_value)],
  late_tmem[, .(Phenotype, variant_id, base_pair_location, p_value)]
))

hs_grn <- merge(hs[CHR == 17], 
                 bel2022_grn[, .(variant_id, base_pair_location)], 
                 by = "variant_id")
late_grn <- merge(late[chr == 17], 
                   bel2022_grn[, .(variant_id, base_pair_location)], 
                   by = "variant_id")
bel2022_grn <- bel2022_grn[variant_id %in% hs_grn$variant_id]
grn <- rbindlist(list(
  bel2022_grn[, .(Phenotype, variant_id, base_pair_location, p_value)],
  hs_grn[, .(Phenotype, variant_id, base_pair_location, p_value)],
  late_grn[, .(Phenotype, variant_id, base_pair_location, p_value)]
))

# colocalization
## TMEM
hs_tmem[, MAF := 1 - (N_case * AF_case + N_ctrl * AF_ctrl) / 
                     (N_case + N_ctrl)]
hs_tmem_coloc <- list(
  MAF=hs_tmem$MAF,
  type="cc",
  beta=hs_tmem$BETA,
  varbeta=hs_tmem$SE^2,
  pvalues=hs_tmem$p_value,
  snp=hs_tmem$variant_id,
  N = hs_pheno[!is.na(hs), .N],
  s = hs_pheno[!is.na(hs), hs_pheno[hs == 1, .N] / .N][1]
)
check_dataset(hs_tmem_coloc)

bel2022_tmem_coloc <- list(
  MAF=bel2022_tmem$effect_allele_frequency,
  type="cc",
  beta=bel2022_tmem$beta,
  varbeta=bel2022_tmem$standard_error^2,
  pvalues=bel2022_tmem$p_value,
  snp=bel2022_tmem$variant_id,
  N = bel2022_tmem[, n_cases + n_controls][1],
  s = bel2022_tmem[, n_cases / (n_cases + n_controls)][1]
)

check_dataset(bel2022_tmem_coloc)

late_tmem[, varbeta := (beta / qnorm(p_value/2, lower.tail = TRUE))^2]
late_tmem_coloc <- list(
  MAF=late_tmem$MAF,
  type="cc",
  beta=late_tmem$beta,
  varbeta=late_tmem$varbeta,
  pvalues=late_tmem$p_value,
  snp=late_tmem$variant_id,
  N = late_pheno[!is.na(late), .N],
  s = late_pheno[!is.na(late), late_pheno[late > 0, .N] / .N]
)
check_dataset(late_tmem_coloc)

coloc.abf(hs_tmem_coloc, bel2022_tmem_coloc, p12 = 1e-5)
coloc.abf(late_tmem_coloc, bel2022_tmem_coloc, p12 = 1e-5)
coloc.abf(hs_tmem_coloc, late_tmem_coloc, p12 = 1e-5)

## GRN
hs_grn[, MAF := 1 - (N_case * AF_case + N_ctrl * AF_ctrl) / 
          (N_case + N_ctrl)]
hs_grn_coloc <- list(
  MAF=hs_grn$MAF,
  type="cc",
  beta=hs_grn$BETA,
  varbeta=hs_grn$SE^2,
  pvalues=hs_grn$p_value,
  snp=hs_grn$variant_id,
  N = hs_pheno[!is.na(hs), .N],
  s = hs_pheno[!is.na(hs), hs_pheno[hs == 1, .N] / .N][1]
)
check_dataset(hs_grn_coloc)

bel2022_grn_coloc <- list(
  MAF=bel2022_grn$effect_allele_frequency,
  type="cc",
  beta=bel2022_grn$beta,
  varbeta=bel2022_grn$standard_error^2,
  pvalues=bel2022_grn$p_value,
  snp=bel2022_grn$variant_id,
  N = bel2022_grn[, n_cases + n_controls][1],
  s = bel2022_grn[, n_cases / (n_cases + n_controls)][1]
)

check_dataset(bel2022_grn_coloc)

late_grn[, varbeta := (beta / qnorm(p_value/2, lower.tail = TRUE))^2]
late_grn_coloc <- list(
  MAF=late_grn$MAF,
  type="cc",
  beta=late_grn$beta,
  varbeta=late_grn$varbeta,
  pvalues=late_grn$p_value,
  snp=late_grn$variant_id,
  N = late_pheno[!is.na(late), .N],
  s = late_pheno[!is.na(late), late_pheno[late > 0, .N] / .N]
)
check_dataset(late_grn_coloc)

coloc.abf(hs_grn_coloc, bel2022_grn_coloc, p12 = 1e-5)
coloc.abf(late_grn_coloc, bel2022_grn_coloc, p12 = 1e-5)
coloc.abf(hs_grn_coloc, late_grn_coloc, p12 = 1e-5)
