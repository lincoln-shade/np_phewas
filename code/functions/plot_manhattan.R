#=============================
# create manhattan plot
#=============================

library(pacman)
p_load(data.table, magrittr, ggplot2, ggrepel)

plot_manhattan <- function(results, signif_only=TRUE, annotate=FALSE) {
  results <- as.data.table(results)
  results[, BP := as.numeric(BP)]
  if (signif_only) {
    results <- results[P < 0.05]
  }
  
  if (length(annotate) > 0) {
    results[SNP %in% annotate, annotate := TRUE]
  }
  
  # cumulative bp length for each chromosome
  chr_bp <- results[, max(BP), CHR] %>%
    setnames("V1", "chr_length") %>%
    # add 20mil buffer to help with spacing of axis labels
    .[, cum_length := cumsum(chr_length) - chr_length + 20000000*(CHR - 1)]
  
  results <- merge(results, chr_bp[, -c("chr_length")], "CHR") %>%
    setorder(CHR, BP) %>%
    .[, bp_cum := BP + cum_length]
  
  x_axis <- results[, (max(bp_cum) + min(bp_cum)) / 2, CHR]
  
  manplot <- results[, ggplot(.SD, aes(bp_cum, -log10(P)))] +
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("#606060", "#000066"), 22)) +
    geom_hline(yintercept = 7.3, color = "black") + #genome-wide significance threshold
    geom_hline(yintercept = 5, color = "black", linetype="dashed") + #suggestibility threshold
    
    scale_x_continuous(label = x_axis$CHR, breaks= x_axis$V1) +
    # remove space between y-axis and baseline
    scale_y_continuous(expand = c(0.05, 0)) +
    
    xlab("Chromosome") +
    ylab("-log(P)") +
    # Add highlighted points
    #geom_point(data=subset(snp_mod, is_highlight=="yes"), color="orange", size=2) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=results[annotate == TRUE], aes(label=SNP), size=5) +
    
    # Custom theme:
    theme_bw() +
    theme(
      # text = element_text(family = "Calibri", size = 20),
      axis.text = element_text(size = 12),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  manplot
}
