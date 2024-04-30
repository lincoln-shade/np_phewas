library(data.table)
library(ggplot2)
library(scales)
dt = np_cpg_rna

# tmem106B expression vs. LATE

p_late_tmem = dt[
  !is.na(TMEM106B) & !is.na(tdp_st4), 
  ggplot(
    .SD, 
    aes(x = factor(tdp_st4), 
        y = log2(TMEM106B) + 1,
        fill = factor(tdp_st4),
        alpha = 0.2))] +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.2) +
  labs(
    x = "LATE pathology stage", 
    y = expression("log"[2] * "(" * italic("TMEM106B") * " expression + 1)")) +
  annotate("text", x = Inf, y = Inf, label = expression(italic("P") * "-value = 0.043"), hjust = 1.1, vjust = 1.1) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = "doc/late_tmem_plot.svg", 
       plot = p_late_tmem, 
       device = "svg", 
       dpi = "retina",
       width = 6,
       height = 3.5,
       units = "in")

p_late_cg09 = dt[
  !is.na(cg09613507) & !is.na(tdp_st4), 
  ggplot(
    .SD, 
    aes(x = factor(tdp_st4), 
        y = cg09613507,
        fill = factor(tdp_st4),
        alpha = 0.2))] +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.2) +
  labs(
    x = "LATE pathology stage", 
    y = expression("cg09613507 % methylated")) +
  annotate("text", x = Inf, y = Inf, label = expression(italic("P") * "-value = 0.0093"), hjust = 1.1, vjust = 1.1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(labels = label_percent(scale = 100))

ggsave(filename = "doc/late_cg09_plot.svg", 
       plot = p_late_cg09, 
       device = "svg", 
       dpi = "retina",
       width = 6,
       height = 3.5,
       units = "in")

p_late_cg23 = dt[
  !is.na(cg23422036) & !is.na(tdp_st4), 
  ggplot(
    .SD, 
    aes(x = factor(tdp_st4), 
        y = cg23422036,
        fill = factor(tdp_st4),
        alpha = 0.2))] +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.2) +
  labs(
    x = "LATE pathology stage", 
    y = expression("cg23422036 % methylated")) +
  annotate("text", x = Inf, y = Inf, label = expression(italic("P") * "-value = 0.10"), hjust = 1.1, vjust = 1.1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(labels = label_percent(scale = 100))

ggsave(filename = "doc/late_cg23_plot.svg", 
       plot = p_late_cg23, 
       device = "svg", 
       dpi = "retina",
       width = 6,
       height = 3.5,
       units = "in")
