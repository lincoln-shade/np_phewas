##--------------------------------------------------
## Single-CpG analysis for genes in gene list
#---------------------------------------------------

library(data.table)
library(synapser)

cpg.meta <- fread("/data_global/ROSMAP/DNA_methylation/ROSMAP_arrayMethylation_metaData.tsv")

genes <- fread("/data_global/UCSC_genome_browser/genes/genes.chr1-22.txt")

gene.list <- fread("02_analysis/gene.list.txt", header = F)
setnames(gene.list, "V1", "Gene")

# text file with one gene name per row
gene.list <- merge(gene.list, genes, by="Gene")

buffer <- 1e5
TargetID <- character()
Gene <- character()
for (i in 1:nrow(gene.list)) {
  gene.name <- gene.list[i, Gene]
  gene.chr <- gene.list[i, Chromosome]
  gene.start <- gene.list[i, Start]
  gene.end <- gene.list[i, End]
  gene.cpg <- cpg.meta[CHR == gene.chr & 
                         MAPINFO > (gene.start - buffer) &
                         MAPINFO < (gene.end + buffer), 
                       TargetID]
  TargetID <- c(TargetID, gene.cpg)
  Gene <- c(Gene, rep(gene.name, length(gene.cpg)))
}
