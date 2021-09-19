##-------------------------
## ordinal regression
##-------------------------

library(pacman)
p_load(data.table, magrittr, MASS, brant)


#.raw file from regression.sh output with variant minor allele values for each participant
raw <- fread("02_analysis/ordinal/top.snps.raw") 

# file with other variables for each participant
load("01_data/data/nacc.ordinal.RData")

# merge to create model matrix file
regression.data <- merge(nacc.ordinal, raw, by = c("FID", "IID"))
n.cols <- ncol(regression.data)

# skip number of columns in nacc.ordinal and the first 4 non-merge-by columns in raw to get to genetic variant columns
skip.col <- ncol(nacc.ordinal) + 4 

# replace colons with periods in SNP names, as colons in variable names messes up polr regression
subColons <- function(x) { 
  x <- gsub(":", ".", x)
}

# add the character "v" in the front of all variant names that do not start with "rs"
# (variant names that start with numbers mess up the regression for some reason)
#logistic.results[-grep("rs", snp), snp := paste0("v", snp)]

for (i in which(colnames(regression.data) %in% 
                colnames(regression.data)[
                  -grep("rs", colnames(regression.data))
                  ][
                    -c(1:skip.col)])
) {
  colnames(regression.data)[i] <- paste0("v", colnames(regression.data)[i])
}

colnames(regression.data) <- subColons(colnames(regression.data))

# response in ordered logistic regression must be a factor
regression.data$NACCARTE <- as.ordered(regression.data$NACCARTE) 

# initialize output vectors
snp <- character(length = (n.cols - skip.col))
brant.p <- numeric(length = (n.cols - skip.col))  
# ordinal regression loop
for (i in (skip.col + 1):n.cols) {
  f = as.formula(paste("NACCARTE ~", colnames(regression.data)[i],
                       " + NPSEX + NACCDAGE + PC1 + PC2 + PC3 + PC4 + PC5 + ",
                       "ADGC.ADC2 + ADGC.ADC3 + ADGC.ADC4 + ADGC.ADC5 + ADGC.ADC6 + ADGC.ADC7", sep = ""))
  m <- polr(f, data = regression.data, Hess = T)
  b <- brant(m)
  snp[(i - skip.col)] <- colnames(regression.data)[i]
  brant.p[i - skip.col] <- b[2, 3]
}

# write outputs to file
brant.test <- data.table(snp, brant.p)
fwrite(brant.test, file = paste0("02_analysis/ordinal/brant.txt"), 
       row.names = F, col.names = T, quote = F, sep = " ")

rm(list = ls())
p_unload(all)
