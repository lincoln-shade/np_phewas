#===========================================================================
# select only Chromosome 19 SNP Matrices for APOE-adjusted analysis
#===========================================================================

# nacc
snpmat_list <- readRDS("data/adc/SNPMatrix/SnpMatrix_list.Rds")
snpmat_list <- snpmat_list[187:192]
saveRDS(snpmat_list, file = "data/adc/SNPMatrix/SnpMatrix_list_chr19.Rds")

# rosmap
snpmat_list <- readRDS("data/rosmap/SNPMatrix/SnpMatrix_list.Rds")
snpmat_list <- snpmat_list[143:146]
saveRDS(snpmat_list, file = "data/rosmap/SNPMatrix/SnpMatrix_list_chr19.Rds")

# ACT
snpmat_list <- readRDS("data/act/SNPMatrix/SnpMatrix_list.Rds")
snpmat_list <- snpmat_list[190:194]
saveRDS(snpmat_list, file = "data/act/SNPMatrix/SnpMatrix_list_chr19.Rds")
