

#======================================================
# annote a variant to the nearest gene
#======================================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-v", "--variant_file", 
                    help = "path to file with variant positions")
parser$add_argument("-g", "--gene_file",  
                    help = "path to file with gene positions")

args <- parser$parse_args(
  c('-v', 'output/gwas/mega/braak56.assoc.logistic',
    '-g', 'raw_data/NCBI38.gene.loc'))

genes <- fread(args$gene_file)
gene_pos <- melt(genes, id.vars = c("V2", "V6"), measure.vars = c("V3", "V4"))
snps <- fread(args$variant_file)

snps_top <- snps[P %in% snps[, min(P), CHR]$V1]

find_closest_gene <- function(chr, bp, gene_pos) {
  chr <- as.character(chr)
  closest_gene <-
    gene_pos[V2 == chr
    ][
      abs(value - bp) == min(abs(value - bp)), V6]
  return(closest_gene)
}

find_closest_gene <- Vectorize(find_closest_gene, 
                               vectorize.args = c("chr", "bp"))

snps_top[, Gene := find_closest_gene(CHR, BP)]
