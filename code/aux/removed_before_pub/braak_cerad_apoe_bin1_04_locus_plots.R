library(data.table)
library(ggplot2)

load("data/mega/coloc/apoe_bin1_coloc_sumstats.Rdata")
braak56_apoe <- coloc_sumstats[[1]]
cerad3_apoe <- coloc_sumstats[[2]]
braak56_bin1 <- coloc_sumstats[[3]]
cerad3_bin1 <- coloc_sumstats[[4]]

braak56_apoe[, Phenotype := 'Braak Stage 5 or 6']
cerad3_apoe[, Phenotype := 'CERAD Score of 3']
braak56_bin1[, Phenotype := 'Braak Stage 5 or 6']
cerad3_bin1[, Phenotype := 'CERAD Score of 3']
apoe <- rbind(braak56_apoe, cerad3_apoe)
bin1 <- rbind(braak56_bin1, cerad3_bin1)

nice_plot <- function(p) {
  p +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(5e-8))
}

p1 <- apoe[, ggplot(.SD, aes(BP, -log10(P), color = Phenotype))]
p2 <- bin1[, ggplot(.SD, aes(BP, -log10(P), color = Phenotype))]


nice_plot(p1) +
  ggtitle("APOE")

nice_plot(p2) +
  ggtitle("BIN1")
