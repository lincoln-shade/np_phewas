##-------------------------
## ordinal regression
##-------------------------
library(pacman)
p_load(data.table, magrittr, rms, stringi)

#-----------------------------------------------
# command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# args[1] = path to text file with covars
# args[2] = path to text file with ordinal phenotypes
# args[3] = path to text file with related ids to remove
# args[4] = phenotype name
# args[5] = snps list file number
#-----------------------------------------------


# .raw file from regression.sh output with variant minor allele values for each participant
snp_list_index <- args[5]
raw <-fread(paste0("data/tmp/ordinal_snp_list_", snp_list_index, ".raw"), header = T)
covar_data <- fread(args[1])
phenotype_data <- fread(args[2], na.strings = "-1")
related_rm <- args[3]
phenotype <- args[4]

# # test args
# raw <- fread(paste0("data/tmp/ordinal_snp_list_", 1, ".raw"), header = T) %>% .[, 1:100]
# covar_data <- fread("data/plink/adc_np.covar")
# phenotype_data <- fread("data/adc_np_ord.txt", na.strings = "-1")
# related_rm <- "data/related_rm/NACCBRAA.remove"
# phenotype <- "NACCBRAA"

# make sure FID and IID are character vectors for merging
covar_data[, `:=`(FID = as.character(FID), IID = as.character(IID))]
raw[, `:=`(FID = as.character(FID), IID = as.character(IID))]
phenotype_data[, `:=`(FID = as.character(FID), IID = as.character(IID))]
set(raw, 
    j = c("PAT", "MAT", "SEX", "PHENOTYPE"), 
    value = list(NULL, NULL, NULL, NULL)
)

# merge to create model matrix file
ord_data <- merge(phenotype_data[, c("FID", "IID", ..phenotype)], 
                  covar_data,
                  by = c("FID", "IID")
)
ord_data <- merge(ord_data, raw, by = c("FID", "IID"))
ord_data[, (phenotype) := as.ordered(get(..phenotype))]

# skip number of columns in covar_data + 1 (the phenotype) get to genetic variant columns
n_cols <- ncol(ord_data)
skip_cols <- ncol(covar_data) + 1L 
  
# replace colons with periods in SNP names, as colons in variable names messes up polr regression
sub_colons <- function(x) { 
  x <- gsub(":", ".", x)
}

# add the character "v" in the front of all variant names that do not start with "rs"
for (i in which(colnames(ord_data) %in% 
                colnames(ord_data)[-grep("rs", colnames(ord_data))
                                   ][-c(1:skip_cols)
                                     ])) {
  colnames(ord_data)[i] <- paste0("v", colnames(ord_data)[i])
}
  
colnames(ord_data) <- sub_colons(colnames(ord_data))

# remove individuals from related_rm folder file
try_fread <- function(path) {
  related_rm_ids <- try(fread(path, fill = T))
  if (inherits(related_rm_ids, "try-error")) {
    return(NULL)
  } else {
    return(related_rm_ids)
  }
}

related_rm_ids <- try_fread(related_rm)

if (!is.null(related_rm_ids)) {
  ord_data <- ord_data[!(IID %in% related_rm_ids$V2)]
}

# check that members of each adgc cohort have phenotype data and
# remove any cohort variables with no members
ord_data <- ord_data[!(is.na(get(phenotype)))]
adgc_cohorts <- unlist(lapply(2:12, function(x) paste0("adgc_", x)))

for (i in adgc_cohorts) {
  if (ord_data[, sum(get(..i))] %in% c(0, ord_data[, .N])) {
    ord_data[, (i) := NULL]
    skip_cols <- skip_cols - 1L
    n_cols <- n_cols - 1L
  }
}

# initialize output table
results <- data.table(
  snp = colnames(ord_data)[(skip_cols+1):n_cols],
  slope = as.numeric(rep(NA, (n_cols - skip_cols))),
  std.er = as.numeric(rep(NA, (n_cols - skip_cols)))
)

# ordinal regression 
run_orm <- function(i) {
  f <- as.formula(paste0(phenotype, " ~ ", colnames(ord_data)[i],
                         " + ", paste(colnames(ord_data)[4:skip_cols], 
                                      collapse = " + "
                                )
                  )
  )
  m <- orm(f, data = ord_data)
  if (m$fail == FALSE) {
    set(
      results,
      i     = i - skip_cols,
      j     = c("slope", "std.er"),
      value = list(
        m$coefficients[colnames(ord_data)[i]],
        m$var[colnames(ord_data)[i], colnames(ord_data)[i]]
      )
    )
  }
}

for (i in (skip_cols + 1L):n_cols) {
  run_orm(i)
}

results[, `:=`(slope = round(slope, 4),
               std.er = round(sqrt(std.er), 4)
          )
]

# write outputs to file
fwrite(results, 
       file = paste0(
         "data/tmp/", 
         phenotype, 
         "_ordinal_results_snp_list_", 
         args[5], 
         ".txt"
       ), 
       quote = F, 
       sep = " "
)
