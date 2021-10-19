library(data.table)
library(ggplot2)
gwas_qtl <- fread("data/tmp/chr2_ENSG00000071575.11_Brain_Cortex_gwas_qtlNPOLD1_bin.tmp")
setnames(gwas_qtl, c("P_gwas", "P_qtl"), c("Cerebral cortex microinfarcts", "TRIB2 expression in cerebral cortex"))
gwas_qtl_long <- melt(gwas_qtl, measure.vars = c("Cerebral cortex microinfarcts", 
                                                 "TRIB2 expression in cerebral cortex"))
setnames(gwas_qtl_long, "variable", "Study")
gwas_qtl_long[, ggplot(.SD, aes(BP_hg38, -log10(value), color = Study))] + 
  geom_point() + 
  theme_minimal() +
  xlab("Chromosome 8 Position") + 
  ylab("-log(P)") + 
  theme(legend.position = c(0.25, 0.8)) +
  theme(legend.text = element_text(size = 10)) + 
  scale_color_manual(values=c("#606060", "#000066"))
