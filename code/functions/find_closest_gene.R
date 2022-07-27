
#======================================================
# annote a variant to the nearest gene
#======================================================

find_closest_gene <- function(chr, bp, g) {
  # chr = chromosome of variant
  # bp = base pair of variant
  # g = data.table with following columns
  #   V2: Chromosome of gene
  #   value: Start or end site of gene (best to use data table with
  #     two rows per gene -- 1 with start and one with end position)
  #   V6: Gene name
  require(data.table)
  g <- as.data.table(g)
  chr <- as.character(chr)
  closest_gene <-
    g[V2 == chr
    ][
      abs(value - bp) == min(abs(value - bp)), V6]
  return(closest_gene)
}

find_closest_gene <- Vectorize(find_closest_gene, 
                               vectorize.args = c("chr", "bp"))