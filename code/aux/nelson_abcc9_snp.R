library(data.table)
library(MASS)

pheno <- fread("data/adc_np_ord.txt")
raw <- fread("data/tmp/abcc9_snp.raw")
covar <- fread("data/plink/adc_np.covar")
data <- merge(pheno[, .(FID, IID, NACCARTE)], raw[, .(FID, IID, `rs1914361:A:G_G`)])
data <- merge(data, covar)

# additive models
none_v_mod <- data[NACCARTE %in% c(0, 2)]
none_v_mod[, NACCARTE := factor(NACCARTE)]
m1 <- none_v_mod[, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]


none_mild_v_mod <- data[NACCARTE %in% c(0, 1, 2)]
none_mild_v_mod[NACCARTE == 1, NACCARTE := 0]
none_mild_v_mod[NACCARTE == 2, NACCARTE := 1]
m2 <- none_mild_v_mod[, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]
summary(m2)

# recessive models
none_v_mod[`rs1914361:A:G_G` == 1, `rs1914361:A:G_G` := 0]
m3 <- none_v_mod[, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]

none_mild_v_mod[`rs1914361:A:G_G` == 1, `rs1914361:A:G_G` := 0]
m4 <- none_mild_v_mod[, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]

# dominance models
none_v_mod <- data[NACCARTE %in% c(0, 2)]
none_v_mod[, NACCARTE := factor(NACCARTE)]
none_v_mod[`rs1914361:A:G_G` == 2, `rs1914361:A:G_G` := 1]
m5 <- none_v_mod[, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]

none_mild_v_mod <- data[NACCARTE %in% c(0, 1, 2)]
none_mild_v_mod[NACCARTE == 1, NACCARTE := 0]
none_mild_v_mod[NACCARTE == 2, NACCARTE := 1]
none_mild_v_mod[`rs1914361:A:G_G` == 2, `rs1914361:A:G_G` := 1]
m6 <- none_mild_v_mod[, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]

# aged over 80 at death
# additive models
none_v_mod <- data[NACCARTE %in% c(0, 2)]
none_v_mod[, NACCARTE := factor(NACCARTE)]
m7 <- none_v_mod[NACCDAGE >= 80, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]


none_mild_v_mod <- data[NACCARTE %in% c(0, 1, 2)]
none_mild_v_mod[NACCARTE == 1, NACCARTE := 0]
none_mild_v_mod[NACCARTE == 2, NACCARTE := 1]
m8 <- none_mild_v_mod[NACCDAGE >= 80, -c("FID", "IID")][, glm(NACCARTE ~ ., data = .SD, family = binomial)]


