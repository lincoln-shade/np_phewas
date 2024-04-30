library(data.table)
library(ggplot2)
library(ordinal)

dt <- fread("data/rosmap_data.csv")
dt[, caa_4gp := ordered(caa_4gp)]
dt[, apoe := factor(apoe, levels = c("apoe33", "apoe22", "apoe23", "apoe34", "apoe44"))]

apoc2_plot = dt[!is.na(rs7247551_G) & !is.na(APOC2), 
                ggplot(.SD, aes(factor(rs7247551_G), sqrt(APOC2)))
]

apoc2_plot + 
  geom_boxplot() +
  labs(
    x = "rs7247551 Genotype",
    y = "sqrt(APOC2 expression)"
  ) +
  scale_x_discrete(labels = c("AA", "AG", "GG")) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust=1, size = 11),
        strip.text.x = element_text(size = 14))

apoe_plot = dt[!is.na(rs7247551_G) & !is.na(APOE), 
                ggplot(.SD, aes(factor(rs7247551_G), log10(APOE)))
]

apoe_plot + 
  geom_boxplot() +
  labs(
    x = "rs7247551 Genotype",
    y = "log10(APOE expression)"
  ) +
  scale_x_discrete(labels = c("AA", "AG", "GG")) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust=1, size = 11),
        strip.text.x = element_text(size = 14))

apoc2_caa_plot = dt[!is.na(caa_4gp) & !is.na(APOC2) & sqrt(APOC2) < 5, 
                    ggplot(.SD, aes(factor(caa_4gp), sqrt(APOC2)))
]

apoc2_caa_plot + 
  geom_boxplot() +
  labs(
    x = "CAA pathology severity",
    y = "sqrt(APOC2 expression)"
  ) +
  scale_x_discrete(labels = c("None", "Mild", "Moderate", "Severe")) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust=1, size = 11),
        strip.text.x = element_text(size = 14))
                    

# regression models: SNP -> gene expression
m_apoe_snp <- dt[, lm(log(APOE) ~ rs7247551_G + msex + age_death + apoe + RIN + pmi)]
m_apoc2_snp <- dt[age_death <80, lm(sqrt(APOC2) ~ rs7247551_G*apoe + msex + age_death + RIN + pmi)]

# gene expression <-> CAA
m_caa_apoe <- dt[, clm(caa_4gp ~ log10(APOE) + msex + age_death + apoe + RIN + pmi)]
m_caa_apoc2 <- dt[, clm(caa_4gp ~ sqrt(APOC2) + msex + age_death + RIN + pmi)]

# SNP -> cpg
m_cg09_snp = dt[, lm(cg09555818 ~ rs7247551_G + msex + age_death + apoe + pmi + conversion.efficiency)]

# cpg <-> gene expression
m_cg09_apoc2 = dt[, lm(sqrt(APOC2) ~ cg09555818 + msex + age_death + apoe + pmi + conversion.efficiency + RIN)]
summary(m_cg09_apoc2)

m_cg13_apoc2 = dt[, lm(sqrt(APOC2) ~ cg13119609 + msex + age_death + apoe + pmi + conversion.efficiency + RIN)]
summary(m_cg13_apoc2)
