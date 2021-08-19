#===============================================================
# Single-variant colocalization analysis
# arg1 should be the GTEx tissue (as it appears in files names)
#===============================================================

library(pacman)
p_load(data.table, magrittr, coloc)
# new
# cargs[1] = p12
# cargs[2] = phenotype_id
# cargs[3] = tissue
# cargsp[4] = chr

cargs <- commandArgs(trailingOnly = T)
p12_cargs <- cargs[1]
out_folder_cargs <- cargs[2]
phenotype_id_cargs <- cargs[3]
tissue_cargs <- cargs[4]
chr_cargs <- cargs[5]

gwas_qtl <- fread(paste0("data/tmp/chr", chr_cargs, "_", phenotype_id_cargs, "_", tissue_cargs, "_gwas_qtl.txt"))
gwas_qtl <- gwas_qtl[!duplicated(gwas_qtl$SNP)]

lead.pos <- gwas_qtl[P_gwas == min(P_gwas), BP_hg38]
lead_snp <- gwas_qtl[BP_hg38 == lead.pos, SNP]
window.radius <- 2e5

gwas_qtl <- gwas_qtl[BP_hg38 > lead.pos - window.radius & BP_hg38 < lead.pos + window.radius]

results <- coloc.abf(
  dataset1 = list(
    beta = gwas_qtl$beta_gwas,
    varbeta = gwas_qtl$var.beta_gwas,
    type = gwas_qtl$type_gwas[1],
    s = gwas_qtl$s[1], 
    N = gwas_qtl$n_gwas, 
    MAF = gwas_qtl$maf_gwas, 
    snp = gwas_qtl$SNP
  ), 
  dataset2 = list(
    beta = gwas_qtl$beta_qtl,
    varbeta = gwas_qtl$var.beta_qtl,
    type = gwas_qtl$type_qtl[1],
    N = gwas_qtl$n_qtl[1],
    MAF = gwas_qtl$maf_qtl,
    snp = gwas_qtl$SNP
  ), 
  p12 = as.numeric(p12_cargs)
)

posteriors <- unname(results$summary)

write(c(cargs[3:5], lead_snp, posteriors), 
      file = paste0(out_folder_cargs, "coloc_results.txt"), 
      ncolumns = 10, 
      append = TRUE)
