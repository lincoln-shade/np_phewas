# prep POLMM results for METAL meta-analysis

library(data.table)
library(argparse)
source("code/functions/strip_alleles.R")
source("code/functions/make_or_95_ci.R")

parser = ArgumentParser()
parser$add_argument("-r", "--results")
parser$add_argument("--bim")
parser$add_argument("-n", "--null_model")
parser$add_argument("--rm_alleles_rsid", default = "FALSE")
parser$add_argument("-o", "--out", default = "./results.csv")
args = parser$parse_args()

append_alleles = function(
    results, 
    bim, 
    snp_col = "SNPID",
    chr_col = "chr",
    a1_col = "A1", 
    a2_col = "A2"
) {
  results_alleles = merge(
    results,
    bim[, c(..chr_col, ..snp_col, ..a1_col, ..a2_col)],
    by = c(chr_col, snp_col)
  )
  
  return(results_alleles)
}

results = fread(args$results)

bim_colnames = c("chr", "SNPID", "CM", "Pos", "A1", "A2")
bim = fread(args$bim, col.names = bim_colnames)

results_al = append_alleles(results, bim)

# add sample sizes to N column
null_model = readRDS(args$null_model)

results_al[, N := null_model$N]

# add standard error estimate column
results_al[, SE := estimate_se_from_beta_and_p(beta, pval.spa)]

# remove alleles from rsID labels if needed
if (as.logical(args$rm_alleles_rsid)) {
  results_al[grep("rs", SNPID), SNPID := strip_alleles(SNPID)]
}


fwrite(results_al, file = args$out)
