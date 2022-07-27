library(ggplot2)
library(data.table)
library(scales)

# GWAS
dtl <- readRDS("output/aux/ord_pwr/ord_pwr_calc1.Rds")
group_colors <- c(`1.4` = "#F8766D", 
                  `1.5` = "#7CAE00", 
                  `1.6` = "#00BFC4", 
                  `1.7` = "#C77CFF")

pw_plot <- dtl[, ggplot(.SD, aes(MAF, Power, color=factor(OR)))] + 
  # geom_smooth(se=FALSE, size = 1.5, span = 1) + 
  geom_line(size = 2) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        legend.justification = c(1, 0),
        legend.position = c(0.9, 0.1),) +
  xlab("Minor Allele Frequency") +
  ylab("Power") +
  labs(color = "Odds Ratio") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = group_colors)

pw_plot

# 4k variants
dtl <- readRDS("output/aux/ord_pwr/ord_pwr_4k_vars.Rds")
group_colors <- c(`1.3` = "#00B6EB", `1.35` = "#A58AFF", `1.4` = "#F8766D", 
                  `1.45` = "#FB61D7", `1.5` = "#7CAE00")

pw_plot <- dtl[, ggplot(.SD, aes(MAF, Power, color=factor(OR)))] + 
  # geom_smooth(se=FALSE, size = 1.5, span = 1) + 
  geom_line(size = 2) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        legend.justification = c(1, 0),
        legend.position = c(0.9, 0.1),) +
  xlab("Minor Allele Frequency") +
  ylab("Power") +
  labs(color = "Odds Ratio") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = group_colors)

pw_plot

# 80 genes
dtl <- readRDS("output/aux/ord_pwr/ord_pwr_80_genes.Rds")
group_colors <- c(`1.2` = "#9590FF", `1.25` = "#E76BF3", `1.3` = "#00B6EB", 
                  `1.35` = "#A58AFF", `1.5` = "#7CAE00")

pw_plot <- dtl[, ggplot(.SD, aes(MAF, Power, color=factor(OR)))] + 
  # geom_smooth(se=FALSE, size = 1.5, span = 1) + 
  geom_line(size = 2) +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        legend.justification = c(1, 0),
        legend.position = c(0.9, 0.1),) +
  xlab("Minor Allele Frequency") +
  ylab("Power") +
  labs(color = "Odds Ratio") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = group_colors)

pw_plot
