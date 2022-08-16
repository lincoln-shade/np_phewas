#===========================================================================
# select only Chromosome 19 SNP Matrices (95-97) for APOE-adjusted analysis
#===========================================================================

snpmat_list <- readRDS("data/mega/SnpMatrix/SnpMatrix_list.Rds")
snpmat_list <- snpmat_list[95:97]
saveRDS(snpmat_list, file = "data/mega/SnpMatrix/SnpMatrix_list95-97.Rds")
