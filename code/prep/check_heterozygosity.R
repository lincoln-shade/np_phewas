#! /usr/bin/Rscript --vanilla
#-----------------------------------------------------------------------------
# Check heterozygosty and create a list of IDs that are > +/- 3 SD from mean
#-----------------------------------------------------------------------------

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-f", "--file",
                    help="PLINK .het file path")
parser$add_argument("-o", "--out",
                    help="Output file path")
args <- parser$parse_args()

het <- fread(args$file)

# (total number genotypes - number homozygous) / total number genotypes
het[, het.rate := (`N(NM)` - `O(HOM)`) / `N(NM)`]

# 3SD threshold for failure
het.fail <- het[het.rate < mean(het.rate) - 3 * sd(het.rate) |
                 het.rate > mean(het.rate) + 3 * sd(het.rate)]
het.fail[, het.dist := (het.rate - mean(het$het.rate)) / sd(het$het.rate)]

write.table(het.fail, file = args$out, quote = F, col.names = F, row.names = F)
