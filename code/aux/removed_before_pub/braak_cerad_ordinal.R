##-------------------------
## ordinal regression
##-------------------------
# use latest version of ordinal package (https://github.com/runehaubo/ordinal)
library(data.table)
library(magrittr)
library(stringi)
library(ordinal)

dt <- readRDS("data/mega/conditional/dt.Rds")
for (i in 1:10) {
  dt[, (paste0("pc", i)) := scale(get(paste0("pc", i)))]
}
dt[, age_death := scale(age_death)]
dt[, age_death2 := scale(age_death2)]
dt[, rs6733839_T := NULL]
dt[, rs4420638_G := NULL]
dt[, `:=`(c('lewy', 'athero', 'grossinf', 'microinf', 'hs', 'braak56', 'cerad3', 'caa', 'arteriol23'),
          c(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL))]
load("data/mega/conditional/covars.Rdata")
covars <- "msex + adgc_2 + adgc_3 + adgc_4 + adgc_5 + adgc_6 + adgc_7 + adgc_8 + adgc_9 + adgc_10 + adgc_11 + adgc_12 + rosmap + ROS + act + ACT2 + ACT3 + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
raw <-fread("data/mega/conditional/bin1.raw", 
            header = T)


# make sure FID and IID are character vectors for merging
dt[, `:=`(FID = as.character(FID), IID = as.character(IID))]
raw[, `:=`(FID = as.character(FID), IID = as.character(IID))]
set(raw, 
    j = c("PAT", "MAT", "SEX", "PHENOTYPE"), 
    value = list(NULL, NULL, NULL, NULL)
)

ord_data <- merge(dt, raw, by = c("FID", "IID"))

# skip number of columns in covar_data + 1 
# (the phenotype) get to genetic variant columns
n_cols <- ncol(ord_data)
skip_cols <- ncol(dt)

# replace colons with periods in SNP names, as colons in variable names 
# messes up regression
sub_colons <- function(x) { 
  x <- gsub(":", ".", x)
}

# add the character "v" in the front of all variant 
# names that do not start with "rs"
for (i in which(colnames(ord_data) %in% 
                colnames(ord_data)[-grep("rs", colnames(ord_data))
                ][-c(1:skip_cols)])) {
  colnames(ord_data)[i] <- paste0("v", colnames(ord_data)[i])
}

colnames(ord_data) <- sub_colons(colnames(ord_data))

# ordinal regression 
run_null_model <- function(phenotype, covariates="") {
  plus_covariates <- ""
  if (covariates != "") {
    plus_covariates <- paste0(" + ", covariates)
  }
  f <- formula(paste0(phenotype, " ~ ", plus_covariates, " + ", covars))
  m <- clm(f, data = ord_data)
}

model_formula <- function(phenotype, snp, covariates="") {
  plus_covariates <- ""
  if (covariates != "") {
    plus_covariates <- paste0(" + ", covariates)
  }
  f <- formula(paste0(phenotype, " ~ ", snp, plus_covariates, " + ", covars))
}

run_clm <- function(i, phenotype, results, covariates='') {
  snp <- colnames(ord_data)[i]
  f <- model_formula(phenotype, snp, covariates)
  m <- clm(f, data = ord_data)
  if (m$convergence$code == 0) {
    sm <- summary(m)
    set(
      results,
      i     = i - skip_cols,
      j     = c("Beta", "SE", "P"),
      value = list(
        sm$coefficients[snp, "Estimate"],
        sm$coefficients[snp, "Std. Error"],
        sm$coefficients[snp, 'Pr(>|z|)']
      )
    )
  }
}

# braak_ord
braak <- data.table(
  SNP = colnames(ord_data)[(skip_cols+1):n_cols],
  Beta = as.numeric(rep(NA, (n_cols - skip_cols))),
  SE = as.numeric(rep(NA, (n_cols - skip_cols))),
  P = as.numeric(rep(NA, (n_cols - skip_cols)))
)

system.time(for (i in (skip_cols + 1L):n_cols) {
  run_clm(i, 'braak_ord', braak)
})

# cerad_ord
cerad <- data.table(
  SNP = colnames(ord_data)[(skip_cols+1):n_cols],
  Beta = as.numeric(rep(NA, (n_cols - skip_cols))),
  SE = as.numeric(rep(NA, (n_cols - skip_cols))),
  P = as.numeric(rep(NA, (n_cols - skip_cols)))
)

null_model <- run_null_model('cerad_ord')
for (i in (skip_cols + 1L):n_cols) {
  run_clm(i, 'cerad_ord', cerad)
}

# braak_ord | cerad
braak_cond <- data.table(
  SNP = colnames(ord_data)[(skip_cols+1):n_cols],
  Beta = as.numeric(rep(NA, (n_cols - skip_cols))),
  SE = as.numeric(rep(NA, (n_cols - skip_cols))),
  P = as.numeric(rep(NA, (n_cols - skip_cols)))
)

for (i in (skip_cols + 1L):n_cols) {
  run_clm(i, 'braak_ord', braak_cond, 'cerad')
}

# cerad_ord | braak
cerad_cond <- data.table(
  SNP = colnames(ord_data)[(skip_cols+1):n_cols],
  Beta = as.numeric(rep(NA, (n_cols - skip_cols))),
  SE = as.numeric(rep(NA, (n_cols - skip_cols))),
  P = as.numeric(rep(NA, (n_cols - skip_cols)))
)

for (i in (skip_cols + 1L):n_cols) {
  run_clm(i, 'cerad_ord', cerad_cond, 'braak')
}

results <- list(braak = braak,
                cerad = cerad,
                braak_cond = braak_cond,
                cerad_cond = cerad_cond)

save(results, file = "data/mega/conditional/braak_cerad_bin1_results.Rdata")


