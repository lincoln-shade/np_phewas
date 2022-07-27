#! /usr/bin/Rscript --vanilla

#-------------------------
# Check heterozygosty
#-------------------------

library(data.table)

het <- fread("tmp/adc_qc1_het.tmp.het")

# (total number genotypes - number homozygous) / total number genotypes
het[, het.rate := (`OBS_CT` - `O(HOM)`) / `OBS_CT`]

# 3SD threshold for failure
het.fail <- het[het.rate < mean(het.rate) - 3 * sd(het.rate) |
                 het.rate > mean(het.rate) + 3 * sd(het.rate)]
het.fail[, het.dist := (het.rate - mean(het$het.rate)) / sd(het$het.rate)]

write.table(het.fail, file = "tmp/adc_qc1_fail_het.tmp", 
            quote = F, col.names = F, row.names = F)
