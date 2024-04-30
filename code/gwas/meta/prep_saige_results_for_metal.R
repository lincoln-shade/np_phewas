# prep SAIGE results for METAL meta-analysis

library(data.table)
library(argparse)
source("code/functions/strip_alleles.R")
source("code/functions/make_or_95_ci.R")

parser = ArgumentParser()
parser$add_argument("-r", "--results")
parser$add_argument("--rm_alleles_rsid", default = "FALSE")
parser$add_argument("-o", "--out", default = "./results.csv")
args = parser$parse_args()

results = fread(args$results)

# add sample sizes to N column
results[, N := N_case + N_ctrl]

# remove alleles from rsID labels if needed
if (as.logical(args$rm_alleles_rsid)) {
  results[grep("rs", MarkerID), MarkerID := strip_alleles(MarkerID)]
}


fwrite(results, file = args$out)
