#-----------------------------------
# creates list of IDs +-3SD of mean
# heterozygosity in PLINK .het file
#-----------------------------------

library(data.table)
cargs <- commandArgs(trailingOnly = TRUE)

het <- fread(cargs[1])

# (total number genotypes - number homozygous) / total number genotypes
het[, het.rate := (`N(NM)` - `O(HOM)`) / `N(NM)`]

# 3SD threshold for failure
het.fail <- het[het.rate < mean(het.rate) - 3 * sd(het.rate) |
                 het.rate > mean(het.rate) + 3 * sd(het.rate)]
het.fail[, het.dist := (het.rate - mean(het$het.rate)) / sd(het$het.rate)]

write.table(het.fail, file = cargs[2], 
            quote = F, col.names = F, row.names = F)
