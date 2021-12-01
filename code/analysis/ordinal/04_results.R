##-----------------------------------------------
## merge ordinal regression results into one dt
##-----------------------------------------------

library(pacman)
p_load(data.table, magrittr, stringi)

cargs <- commandArgs(trailingOnly = TRUE)
out_prefix <- "output/ordinal"
pheno <- "NPOLD2" # cargs[1]
# include trailing "/"
directory <- "data/tmp/"

lists <- list.files(directory) %>% 
  .[grep(paste0(pheno, "_ordinal_results"), .)]

list_set <- vector("list", length=length(lists))

for (i in 1:length(lists)) {
  list_set[[i]] <- fread(paste0(directory, lists[i]), fill = TRUE)
}

ordinal_results <- rbindlist(list_set)
ordinal_results <- ordinal_results[!(is.na(Beta))]

# create columns for P value and OR for ordinal regression
ordinal_results[, STAT := round(Beta / SE, 4)]
ordinal_results[, P := signif(pnorm(abs(STAT), lower.tail = F) * 2, 4)]
ordinal_results[, OR := round(exp(Beta), 2)]

# rename SNPs in ordinal regression to match those from logistic
ordinal_results[, SNP := stri_replace_last_regex(SNP, "_[:alpha:]*", "")]
ordinal_results[, SNP := stri_replace_all_fixed(SNP, ".", ":")]
ordinal_results[, SNP := stri_replace_first_fixed(SNP, "v", "")]

setorder(ordinal_results, P)
print(head(ordinal_results))

# setnames(results, 
#          c("OR", "STAT", "P", "L95", "U95", "SE"), 
#          c("OR.Log", "STAT.Log", "P.Log", "L95.Log", "U95.Log", "SE.Log"))
# setcolorder(results, 
#             c("CHR", "SNP", "BP", "A1", "OR.Ord", "SE.Ord", "P.Ord", "STAT.Ord", 
#               "OR.Log", "P.Log", "STAT.Log", "L95.Log", "U95.Log"))

if (!(dir.exists(out_prefix))) {
  dir.create(out_prefix)
}

fwrite(ordinal_results, 
       file = paste0(out_prefix, "/", pheno, "_ordinal_results.txt"), 
       quote = F,
       na = "NA",
       sep = " ")
