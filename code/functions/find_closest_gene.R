
#======================================================
# annote a variant to the nearest gene
#======================================================

find_closest_gene <- function(chr, bp, g, 
                              chrom_col="chr", 
                              pos_col="value",
                              gene_col="gene_name") {
  # chr = chromosome of variant
  # bp = base pair of variant
  # g = data.table with following columns
  #   chrom_col: Chromosome of gene
  #   pos_col: Start or end site of gene (best to use data table with
  #     two rows per gene -- 1 with start and one with end position)
  #   gene_col: Gene name
  require(data.table)
  g <- as.data.table(g)
  snp_chr <- as.character(chr)
  closest_gene <-
    g[as.character(get(chrom_col)) == snp_chr][
      abs(get(pos_col) - bp) == min(abs(get(pos_col) - bp)), get(gene_col)]
  
  if (length(closest_gene) > 1) {
    closest_gene <- paste(closest_gene, collapse = "/")
  }
  
  return(closest_gene)
}

find_closest_gene <- Vectorize(find_closest_gene, 
                               vectorize.args = c("chr", "bp"))