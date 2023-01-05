library(data.table)
library(ggplot2)
library(ggthemes)

hs <- fread("output/gwas/mega/saige/hs_saig_results.txt")
late <- fread("output/gwas/mega/polmm/late_polmm_results.txt")
bell2022 <- 
  fread("/data_global/summary_statistics/GCST90027158_buildGRCh38.tsv.gz")
radius <- 2e5L
tmem_lead_snp <- hs[CHR == 7][p.value == min(p.value), MarkerID]
grn_lead_snp <- hs[CHR == 17][p.value == min(p.value), MarkerID]
tmem_lead_pos <- bell2022[chromosome == 7
  ][variant_id == tmem_lead_snp, base_pair_location]
grn_lead_pos <- bell2022[chromosome == 17
  ][variant_id == grn_lead_snp, base_pair_location]

hs[, Phenotype := "HS"]
late[, Phenotype := "LATE-NC"]
bell2022[, Phenotype := "AD"]

bell2022[, p_value := as.numeric(p_value)]

setnames(hs, c("MarkerID", "p.value"), c("variant_id", "p_value"))
setnames(late, c("SNPID", "pval.spa"), c("variant_id", "p_value"))

bell2022_tmem <- 
  bell2022[chromosome == 7
    ][(base_pair_location >= tmem_lead_pos - radius) &
        (base_pair_location <= tmem_lead_pos + radius)]

bell2022_grn <- 
  bell2022[chromosome == 17
  ][(base_pair_location >= grn_lead_pos - radius) &
      (base_pair_location <= grn_lead_pos + radius)]

hs_tmem <- 
  hs[CHR == 7
    ][variant_id %in% bell2022_tmem[, c(variant_id, variant_alternate_id)]]
late_tmem <- 
  late[chr == 7
    ][variant_id %in% bell2022_tmem[, c(variant_id, variant_alternate_id)]]
bell2022_tmem <- 
  bell2022_tmem[variant_id %in% hs_tmem$variant_id
    ][!which(duplicated(variant_id))]

hs_tmem <- merge(hs_tmem, 
                 bell2022_tmem[, .(variant_id, base_pair_location)],
                 by = "variant_id")
late_tmem <- merge(late_tmem, 
                 bell2022_tmem[, .(variant_id, base_pair_location)],
                 by = "variant_id")

tmem <- rbindlist(list(
  bell2022_tmem[, .(Phenotype, base_pair_location, p_value)],
  hs_tmem[, .(Phenotype, base_pair_location, p_value)],
  late_tmem[, .(Phenotype, base_pair_location, p_value)]))

tmem_ggp <- 
  tmem[, 
       ggplot(.SD, 
              aes(base_pair_location / 1e6L, 
                  -log10(p_value), 
                  color = Phenotype))]

tmem_ggp <- 
  tmem_ggp + 
  geom_point() +
  labs(x = "Chromosome 7 Position (MB)",
       y = "-log(P)") +
  theme_minimal() +
  theme(text = element_text(size = 16)) +  
  scale_color_colorblind()

ggsave(plot = tmem_ggp, 
     file = "doc/fig4a_tmem_coloc.png",
     units = "in",
     width = 10,
     height = 8)

# GRN
hs_grn <- 
  hs[CHR == 17
  ][variant_id %in% bell2022_grn[, c(variant_id, variant_alternate_id)]]
late_grn <- 
  late[chr == 17
  ][variant_id %in% bell2022_grn[, c(variant_id, variant_alternate_id)]]
bell2022_grn <- 
  bell2022_grn[variant_id %in% hs_grn$variant_id
  ][!which(duplicated(variant_id))]

hs_grn <- merge(hs_grn, 
                 bell2022_grn[, .(variant_id, base_pair_location)],
                 by = "variant_id")
late_grn <- merge(late_grn, 
                   bell2022_grn[, .(variant_id, base_pair_location)],
                   by = "variant_id")

grn <- rbindlist(list(
  bell2022_grn[, .(Phenotype, base_pair_location, p_value)],
  hs_grn[, .(Phenotype, base_pair_location, p_value)],
  late_grn[, .(Phenotype, base_pair_location, p_value)]))

grn_ggp <- 
  grn[, 
       ggplot(.SD, 
              aes(base_pair_location / 1e6L, 
                  -log10(p_value), 
                  color = Phenotype))]

grn_ggp <- 
  grn_ggp + 
  geom_point() +
  labs(x = "Chromosome 17 Position (MB)",
       y = "-log(P)") +
  theme_minimal() +
  theme(text = element_text(size = 16)) +  
  scale_color_colorblind()

ggsave(plot = grn_ggp, 
       file = "doc/fig4b_grn_coloc.png",
       units = "in",
       width = 10,
       height = 8)
