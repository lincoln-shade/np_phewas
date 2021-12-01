#==================================
# categorize the NACC NP variables
# by outcome type, coding, etc.
#==================================

library(data.table)
library(readxl)
library(magrittr)

act <- as.data.table(read_xlsx("raw_data/act/DataForE235_20210421.xlsx",
                               sheet = 2))

adc <- fread("data/plink/adc_np.pheno")
