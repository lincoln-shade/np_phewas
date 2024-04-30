library(data.table)
library(ordinal)
library(mediation)
library(MASS)

dt <- fread("data/rosmap_data.csv")
dt[, caa_4gp := ordered(caa_4gp)]
dt[, apoe := factor(apoe, levels = c("apoe33", "apoe22", "apoe23", "apoe34", "apoe44"))]

m_null <- dt[, clm(caa_4gp ~ scale(conversion.efficiency) + scale(pmi) + 
                       msex + scale(age_death) + study + apoe)]
m1 <- dt[, clm(caa_4gp ~ scale(cg04401876) + scale(cg09555818) + 
                   scale(cg10169327) + scale(cg13119609) + 
                   scale(conversion.efficiency) + scale(pmi) + msex + 
                   scale(age_death) + study + apoe)]
m2 <- dt[, clm(caa_4gp ~ scale(cg04401876) + scale(conversion.efficiency) 
               + scale(pmi) + msex + scale(age_death) + study + apoe)]
m3 <- dt[, clm(caa_4gp ~ scale(cg09555818) + scale(conversion.efficiency) 
               + scale(pmi) + msex + scale(age_death) + study + apoe)]
m4 <- dt[, clm(caa_4gp ~ scale(cg10169327) + scale(conversion.efficiency) + 
                   scale(pmi) + msex + scale(age_death) + study + apoe)]
m5 <- dt[, clm(caa_4gp ~ scale(cg13119609) + scale(conversion.efficiency) + 
                   scale(pmi) + msex + scale(age_death) + study + apoe)]
anova(m1, m_null)
anova(m2, m_null) # cg04401876

# # apoe vs methylation
# dt[apoe %in% c("apoe34", "apoe44"), apoe4 := 1]
# dt[!is.na(apoe) & is.na(apoe4), apoe4 := 0]
# cpg_m0 <- dt[!is.na(apoe), lm(cg09555818 ~ 1)]
# cpg_m <- dt[, lm(cg04401876 ~ apoe)]

# CAA ~ SNP
m_snp <- dt[, clm(caa_4gp ~ rs7247551_G + scale(conversion.efficiency) + scale(pmi) + msex + 
                      scale(age_death) + study + apoe)]
m1_snp <- dt[, clm(caa_4gp ~  rs7247551_G +
                       scale(cg09555818) + 
                       scale(cg13119609) + 
                       scale(conversion.efficiency) + scale(pmi) + msex + 
                       scale(age_death) + study + apoe)]
m2_snp <- dt[, clm(caa_4gp ~  rs7247551_G +
                       scale(cg09555818) +
                       scale(conversion.efficiency) + scale(pmi) + msex + 
                       scale(age_death) + study + apoe)]
anova(m_snp, m2_snp)

# mRNA expression models
m_apoe_null <- dt[, lm(log(APOE) ~ msex + age_death + apoe + RIN + pmi + conversion.efficiency)]
m_apoe_snp <- dt[, lm(log(APOE) ~ rs7247551_G + msex + age_death + apoe + RIN + pmi + conversion.efficiency)]
m_apoe_snp_cpg <- dt[, lm(log(APOE) ~ rs7247551_G + cg04401876 + cg09555818 + cg10169327 + cg13119609 + conversion.efficiency + msex + age_death + apoe + RIN + pmi)]
m_apoe_pc1 <- dt[, lm(log(APOE) ~ PC1 + conversion.efficiency + msex + age_death + apoe + RIN + pmi)]

m_apoc2_null <- dt[APOC2 < 10 & !is.na(rs7247551_G), lm(sqrt(APOC2) ~ msex + age_death + apoe + RIN + pmi + conversion.efficiency)]
m_apoc2_snp <- dt[APOC2 < 10, lm(sqrt(APOC2) ~ rs7247551_G + msex + age_death + apoe + RIN + pmi + conversion.efficiency)]
m_apoc2_snp_cpg <- dt[APOC2 < 10, lm(sqrt(APOC2) ~ rs7247551_G + cg04401876 + cg09555818 + cg10169327 + cg13119609 + conversion.efficiency + msex + age_death + apoe + RIN + pmi)]
m_apoc2_cpg <- dt[APOC2 < 10 & !is.na(rs7247551_G), lm(sqrt(APOC2) ~ cg04401876 + cg09555818 + cg10169327 + cg13119609 + conversion.efficiency + msex + age_death + apoe + RIN + pmi)]
m_apoc2_pc1 <- dt[APOC2 < 10 & !is.na(rs7247551_G), lm(sqrt(APOC2) ~ PC1 + conversion.efficiency + msex + age_death + apoe + RIN + pmi)]
m_null <- dt[!is.na(rs7247551_G), clm(caa_4gp ~  scale(conversion.efficiency) + scale(pmi) + msex + 
                                         +                       scale(age_death) + study + apoe)]

# mRNA -> CAA models
# APOE
m_caa_apoe_null <- dt[!is.na(rs7247551_G) & !is.na(APOE), clm(caa_4gp ~ msex + age_death + apoe + RIN + pmi)]
m_caa_apoe <- dt[!is.na(APOE), clm(caa_4gp ~ sqrt(APOE) + msex + age_death + apoe + RIN + pmi)]
m_caa_apoe_snp <- dt[!is.na(rs7247551_G) & !is.na(APOE), clm(caa_4gp ~ log(APOE) + rs7247551_G + msex + age_death + apoe + RIN + pmi)]
m_caa_apoe_snp_pc1 <- dt[!is.na(rs7247551_G) & !is.na(APOE), clm(caa_4gp ~ log(APOE) + rs7247551_G + PC1 + msex + age_death + apoe + RIN + pmi)]

# APOC2
m_caa_apoc2_null <- dt[!is.na(rs7247551_G) & !is.na(APOC2), clm(caa_4gp ~ msex + age_death + apoe + RIN + pmi)]
m_caa_apoc2 <- dt[!is.na(rs7247551_G) & !is.na(APOC2), clm(caa_4gp ~ APOC2 + msex + age_death + apoe + RIN + pmi)]
m_caa_apoc2_snp <- dt[!is.na(rs7247551_G) & !is.na(APOC2), clm(caa_4gp ~ sqrt(APOC2) + rs7247551_G + msex + age_death + apoe + RIN + pmi)]
m_caa_apoc2_snp_pc1 <- dt[!is.na(rs7247551_G) & !is.na(APOC2), clm(caa_4gp ~ sqrt(APOC2) + rs7247551_G + PC1 + msex + age_death + apoe + RIN + pmi)]


#plots
dt[!is.na(rs7247551_G), ggplot(.SD, aes(factor(rs7247551_G), sqrt(APOC2))) + geom_boxplot()]
dt[!is.na(rs7247551_G), ggplot(.SD, aes(factor(rs7247551_G), sqrt(APOE))) + geom_boxplot()]
dt[, ggplot(.SD, aes(PC1, sqrt(APOE))) + geom_point()]
dt[sqrt(APOC2) < 2.5, ggplot(.SD, aes(PC1, sqrt(APOC2))) + geom_point()]
