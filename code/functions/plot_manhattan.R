#=============================
# create manhattan plot
#=============================

plot_manhattan <- function(results, 
                           thin=TRUE, 
                           annotate_snps=FALSE,
                           annotate_labels=FALSE,
                           highlight=FALSE,
                           font_size=12,
                           chr_col = "chr",
                           pos_col = "POS",
                           id_col = "SNPID",
                           pval_col = "pval.spa",
                           max_yval = 20) {
  require(data.table)
  require(magrittr)
  require(ggplot2)
  require(ggrepel)
  
  results <- as.data.table(results)
  setnames(results, 
           c(chr_col, pos_col, id_col, pval_col), 
           c('CHR', 'BP', 'SNP', "P"))
  setorder(results, CHR, BP)
  results[, BP := as.numeric(BP)]
  if (thin) {
    results_signif <- results[P < 0.05]
    results_nonsignif <- results[P >= 0.05][sample(.N, .N / 10)]
    results <- rbind(results_signif, results_nonsignif)
    setorder(results, CHR, BP)
  }
  
  results[, annotate := FALSE]
  results[, is_highlight := FALSE]
  results[, Label := ""]
  
  # cumulative bp length for each chromosome
  chr_bp <- results[, max(BP), CHR] %>%
    setnames("V1", "chr_length") %>%
    # add 20mil buffer to help with spacing of axis labels
    .[, cum_length := cumsum(chr_length) - chr_length + 20000000*(CHR - 1)]
  
  results <- merge(results, chr_bp[, -c("chr_length")], "CHR") %>%
    setorder(CHR, BP) %>%
    .[, bp_cum := BP + cum_length]
  
  # add annotation labels and highlights, if applicable
  if (is.character(annotate_snps)) {
    results[SNP %in% annotate_snps, 
            annotate := TRUE]
    if (is.character(annotate_labels) & 
        length(annotate_labels) == length(annotate_snps)) {
      results[annotate == TRUE, Label := annotate_labels]
    } else {
      results[annotate == TRUE, Label := annotate_snps]
    }
    
    if (highlight == TRUE) {
      radius = 5e4L
      for (snp in annotate_snps) {
        chr <- results[SNP == snp, CHR]
        bp <- results[SNP == snp, bp_cum]
        results[CHR == chr & bp_cum > bp - radius & bp_cum < bp + radius,
                is_highlight := TRUE]
      }
    }
  }
  
  results[, neg_logp := -log10(P)]
  results[neg_logp > max_yval, neg_logp := max_yval]
  x_axis <- results[, (max(bp_cum) + min(bp_cum)) / 2, CHR]
  max_y <- results[, min(max(-log10(P)), max_yval)]
  manplot <- results[, ggplot(.SD, aes(bp_cum, neg_logp))] +
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("#606060", "#000066"), 22)) +
    geom_hline(yintercept = 7.3, color = "black") + #genome-wide significance threshold
    geom_hline(yintercept = 5, color = "black", linetype="dashed") + #suggestibility threshold
    
    scale_x_continuous(label = x_axis$CHR, breaks= x_axis$V1) +
    # remove space between y-axis and baseline
    # scale_y_continuous(expand = c(0.05, 0))+ 
    scale_y_continuous(breaks = seq(2, max_y, 2)) +
    
    xlab("Chromosome") +
    ylab("-log(P)") +
    # Add highlighted points
    geom_point(data=results[is_highlight == TRUE], color="orange", size=2) +
    
    # Add label
    geom_text(data=results[annotate == TRUE],
              aes(label=Label), 
              size=4, 
              family = "Arial",
              fontface="italic",
              vjust = -1) +
    
    # Custom theme:
    theme_bw() +
    theme(
      # text = element_text(family = "Calibri", size = 20),
      axis.text = element_text(size = font_size),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  manplot
}
