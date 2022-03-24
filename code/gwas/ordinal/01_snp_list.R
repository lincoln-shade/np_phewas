##---------------------------------------------------------------------
## make lists of snps of length 30K for analysis blocks 
##---------------------------------------------------------------------

library(pacman)
p_load(data.table, magrittr)
# args <- commandArgs(trailingOnly = TRUE)

# load snps from .bim file
bim_file <- "data/plink/adc_np.bim"# args[1]
bim <- fread(bim_file, header = F) %>% 
  .[, .(V1, V2)] %>% 
  setnames(., c("V1", "V2"), c("chr", "snp")) %>% 
  .[, chr := factor(chr)]

raw_length <- 30000 # desired approximate number of variants for each .raw file
n_lists <- nrow(bim) %/% raw_length
bim[, list_n := rep(1:n_lists, length.out=.N)]

for (i in seq_len(n_lists)) {
  fwrite(bim[list_n == i, .(snp)], 
         file = paste0("data/tmp/ordinal_snp_list_", i, ".tmp"), 
         row.names = F, 
         col.names = F, 
         quote = F, 
         sep = " "
  )
}

list_index <- data.table(1:n_lists)
fwrite(list_index, file = "data/tmp/snp_list_index.tmp", row.names = F, col.names = F, quote = F, sep = " ")
