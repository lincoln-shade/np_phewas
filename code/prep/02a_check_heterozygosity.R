#-------------------------
# Check heterozygosty
#-------------------------

library(data.table)

het <- fread("data/tmp/nacc_qc1_het.tmp.het")

# (total number genotypes - number homozygous) / total number genotypes
het[, het.rate := (`N(NM)` - `O(HOM)`) / `N(NM)`]

# 3SD threshold for failure
het.fail <- het[het.rate < mean(het.rate) - 3 * sd(het.rate) |
                 het.rate > mean(het.rate) + 3 * sd(het.rate)]
het.fail[, het.dist := (het.rate - mean(het$het.rate)) / sd(het$het.rate)]

write.table(het.fail, file = "data/tmp/nacc_qc1_fail_het.tmp.txt", 
            quote = F, col.names = F, row.names = F)
