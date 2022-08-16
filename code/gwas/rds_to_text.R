#===================================
# write text files of POLMM results
#===================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-r", "--rds", help = "results .Rds file")
parser$add_argument("-o", "--out", help = "output file name")
args <- parser$parse_args()

results <- readRDS(args$rds)
fwrite(results, file = args$out, sep = " ", quote = FALSE)