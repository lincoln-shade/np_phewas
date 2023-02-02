library(data.table)
library(gwasglue)
library(coloc)

calc_betavar_from_beta_pval <- function(beta, pval) {
  z <- qnorm(pval / 2, lower.tail = FALSE)
  se <- beta / z
  betavar <- se^2
  betavar
}
calc_betavar_from_beta_pval <- Vectorize(calc_betavar_from_beta_pval)

make_coloc_dataset <- function(l) {
  ds <- c(l$dataset, l$metadata)
}

braak_bin1 <- readRDS("output/coloc/mega/braak/GWAS/gwas_chr2_127133851_A_C_b38.Rds")
gtex_bin1 <- readRDS("output/coloc/mega/braak/GTEx/Brain_Cerebellum_ENSG00000136717.14_chr2_127133851_A_C_b38.Rds")
braak_results <- fread("output/gwas/mega/polmm/braak_polmm_results.txt")
braak_results <- braak_results[chr == 2][SNPID %in% braak_bin1$dataset$rsID]
braak_bin1$dataset <- merge(braak_bin1$dataset,
                            braak_results[, .(SNPID, beta)],
                            by.x = "rsID",
                            by.y = "SNPID")
gtex_bin1$dataset <- gtex_bin1$dataset[snp %in% braak_bin1$dataset$chr_bp_ref_alt_b38]
braak_bin1$dataset <- braak_bin1$dataset[chr_bp_ref_alt_b38 %in% gtex_bin1$dataset$snp]
snps <- braak_bin1$dataset$rsID

ld_mat <- ld_matrix_local(snps, 
                with_alleles = FALSE, 
                bfile = "/home/lmsh224/1000g_plink/EUR",
                plink_bin = genetics.binaRies::get_plink_binary())
keep_snps <- colnames(ld_mat)
ld_mat_gwas <- ld_matrix_local(keep_snps, 
                          with_alleles = FALSE, 
                          bfile = "data/mega/mega_np",
                          plink_bin = genetics.binaRies::get_plink_binary())

braak_bin1$dataset <- braak_bin1$dataset[rsID %in% keep_snps]
gtex_bin1$dataset <- gtex_bin1$dataset[snp %in% braak_bin1$dataset$chr_bp_ref_alt_b38]

braak_bin1$metadata$LD <- ld_mat_gwas
gtex_bin1$metadata$LD <- ld_mat

gtex_bin1$dataset <- merge(gtex_bin1$dataset,
                           braak_bin1$dataset[, .(rsID, chr_bp_ref_alt_b38)],
                           by.x = "snp",
                           by.y = "chr_bp_ref_alt_b38")
gtex_bin1$dataset$snp <- gtex_bin1$dataset$rsID

braak_bin1$dataset[, varbeta := calc_betavar_from_beta_pval(beta, pvalues)]
braak_bin1$dataset[, snp := rsID]
braak_bin1$dataset <- braak_bin1$dataset[order(match(braak_bin1$dataset$rsID, colnames(braak_bin1$metadata$LD)))]
gtex_bin1$dataset <- gtex_bin1$dataset[order(match(gtex_bin1$dataset$snp, colnames(gtex_bin1$metadata$LD)))]

d1 <- make_coloc_dataset(braak_bin1)
d2 <- make_coloc_dataset(gtex_bin1)
check_dataset(d1, req = "LD")
check_dataset(d2, req = "LD")

s1 <- runsusie(d1)
s2 <- runsusie(d2)

if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(s1, s2)
  print(susie.res$summary)
}


