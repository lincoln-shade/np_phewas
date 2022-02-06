
library(data.table)
library(magrittr)

adc_pheno <- fread("data/adc/adc_np.pheno", na.strings = '-1')

shared_phenos <- c('HS', 'NACCMICRO', 'ASC', 'NACCINF', 'BRAAK', 'NEUR', 'LEWY')
ibd <- fread("data/rosmap/rosmap_adc.kin0")

rm_ids <- ibd[, .(FID2, ID2)]
fwrite(rm_ids, file = "data/rosmap/rosmap_adc_related.txt",
       quote = FALSE, sep = ' ', col.names = FALSE)
