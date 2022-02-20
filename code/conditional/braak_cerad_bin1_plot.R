library(data.table)
library(magrittr)
library(ggplot2)
library(stringi)
library(scales)
library(ggrepel)

load("data/mega/conditional/braak_cerad_bin1_results.Rdata")
snps <- fread("data/mega/mega_np.bim")
setnames(snps, colnames(snps), c('CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'))

results$braak[, Phenotype := 'Braak Stage']
results$cerad[, Phenotype := 'CERAD Score']
results$braak_cond[, Phenotype := 'Braak Stage']
results$cerad_cond[, Phenotype := 'CERAD Score']
uncond <- rbind(results$braak, results$cerad)
cond <- rbind(results$braak_cond, results$cerad_cond)

rm_A1 <- function(snps) {
  snps <- stri_replace_all_regex(snps, "_[ACTG]*", "")
}

uncond[, SNP := rm_A1(SNP)]
cond[, SNP := rm_A1(SNP)]

uncond <- merge(uncond, snps[CHR == 2], 'SNP')
cond <- merge(cond, snps[CHR == 2], 'SNP')

top_snp <- 'rs6733839'
uncond[SNP %in% top_snp, annot := TRUE]
cond[SNP %in% top_snp, annot := TRUE]

uncond[, Analysis := 'Unconditional']
cond[, Analysis := 'Conditional']

res <- rbind(uncond, cond)
res[, Analysisf := factor(Analysis, levels = c('Unconditional', 'Conditional'))]

nice_plot <- function(p) {
  p +
    geom_point() +
    theme_minimal() +
    geom_hline(yintercept = -log10(5e-8))
}

p1 <- uncond[, ggplot(.SD, aes(BP, -log10(P), color = Phenotype))]
p2 <- cond[, ggplot(.SD, aes(BP, -log10(P), color = Phenotype))]
p3 <- res[, ggplot(.SD, aes(BP, -log10(P), color = Phenotype))]



nice_plot(p3) +
  ggtitle("BIN1 Regional Association Plots of Braak Stage and CERAD Score") +
  scale_y_continuous(limits = c(0, 11), breaks = 1:10) +
  scale_x_continuous(label=comma) +
  geom_label_repel(data = res[annot == TRUE], aes(label = SNP)) + 
  facet_grid(cols = vars(Analysisf)) +
  theme(text = element_text(size = 16),
        legend.justification = c(1, 0),
        legend.position = c(0.5, 0.4),)
