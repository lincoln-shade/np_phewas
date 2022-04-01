#! /usr/bin/Rscript --vanilla

#===============================================================
# Single-variant colocalization analysis
# arg1 should be the GTEx tissue (as it appears in files names)
#===============================================================

library(data.table)
library(magrittr)
library(coloc)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-f", "--file", 
                    help="GWAS and QTL merged sumstats file")
parser$add_argument("-o", "--out", 
                    help="output file path")
parser$add_argument("--p12", help="prior probability of colocalization")
parser$add_argument("-q", "--qtl_pheno", 
                    help="QTL phenotype")
parser$add_argument("-t", "--tissue", 
                    help="QTL tissue")
parser$add_argument("--chr", help="chromosome integer")
args <- parser$parse_args()

# #########################################
# args <- list(
#   file = "tmp/chr20_ENSG00000132669.12_Thyroid_gwas_qtlarteriol23.tmp",
#   out = "tmp/arteriol23_coloc_results.txt",
#   p12 = "0.00001",
#   qtl_pheno = "ENSG00000132669.12",
#   tissue = "Thyroid",
#   chr = "20"
# )
# #########################################

p12 <- as.numeric(args$p12)
gwas_qtl <- fread(args$file)
gwas_qtl <- gwas_qtl[!duplicated(gwas_qtl$SNP)]

lead.pos <- gwas_qtl[P_gwas == min(P_gwas), BP_hg38][1]
lead_snp <- gwas_qtl[BP_hg38 == lead.pos, SNP]
window.radius <- 2e5

gwas_qtl <- gwas_qtl[BP_hg38 > lead.pos - window.radius & 
                       BP_hg38 < lead.pos + window.radius]


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
  p12 = p12
)



# write output
posteriors <- formatC(unname(results$summary))
top_snp <- gwas_qtl[SNP == lead_snp]
cols <- c('Phenotype', 'Tissue', 'Chromosome', 'SNP', 'beta_GWAS', 'P_GWAS',
          'beta_QTL', 'P_QTL', 'NSNP', 'PPH0', 'PPH1', 'PPH2', 'PPH3', 
          'PPH4')
if (!(file.exists(args$out))) {
  write(cols,
        file = args$out,
        ncolumns = length(cols))
}
write(c(args$qtl_pheno, args$tissue, args$chr, lead_snp, top_snp$beta_gwas,
        top_snp$P_gwas, top_snp$beta_qtl, top_snp$P_qtl, posteriors), 
      file = args$out, 
      ncolumns = length(cols), 
      append = TRUE)
