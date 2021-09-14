#===============================================================
# Single-variant colocalization analysis
# arg1 should be the GTEx tissue (as it appears in files names)
#===============================================================

library(pacman)
p_load(data.table, magrittr, coloc)
# new
# args[1] = p12
# args[2] = phenotype_id
# args[3] = tissue
# argsp[4] = chr

args <- commandArgs(trailingOnly = T)
p12_args <- args[1]
out_folder_args <- args[2]
gwas_phenotype <- args[3]
phenotype_id_args <- args[4]
tissue_args <- args[5]
chr_args <- args[6]

gwas_qtl <- fread(paste0("data/tmp/chr", chr_args, "_", phenotype_id_args, "_", tissue_args, "_gwas_qtl", gwas_phenotype, ".tmp"))
gwas_qtl <- gwas_qtl[!duplicated(gwas_qtl$SNP)]

lead.pos <- gwas_qtl[P_gwas == min(P_gwas), BP_hg38]
lead_snp <- gwas_qtl[BP_hg38 == lead.pos, SNP]
window.radius <- 2e5

gwas_qtl <- gwas_qtl[BP_hg38 > lead.pos - window.radius & 
                       BP_hg38 < lead.pos + window.radius
]


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
  p12 = as.numeric(p12_args)
)

posteriors <- unname(results$summary)

write(c(args[4:6], lead_snp, posteriors), 
      file = paste0(out_folder_args, "_coloc_results.txt"), 
      ncolumns = 10, 
      append = TRUE)
