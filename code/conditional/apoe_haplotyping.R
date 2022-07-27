
library(data.table)
library(stringi)

nacc <- fread('data/adc/apoe_snps.raw')
rosmap <- fread('data/rosmap/apoe_snps.raw')
act <- fread('data/act/apoe_snps.raw')
adni <- fread('data/adni/apoe_snps.raw')
covar <- fread('data/mega/mega_np.covar')

nacc[, FID := as.character(FID)]
nacc[, FID := '0']
nacc[, IID := stri_replace_first_regex(IID, '.*_', '')]
setnames(nacc, colnames(nacc)[7:8], c('rs429358_C', 'rs7412_T'))

rosmap[, FID := as.character(FID)]
rosmap[, FID := '1']

act[, FID := as.character(FID)]
act[, FID := '2']
act[, IID := stri_replace_first_regex(IID, '.*_', '')]
setnames(act, colnames(act)[7:8], c('rs429358_C', 'rs7412_T'))

adni[, FID := as.character(FID)]
adni[, FID := 'ADNI']

apoe <- rbindlist(list(nacc, rosmap, act, adni), use.names = TRUE)
apoe[, `:=`(PAT = NULL, MAT = NULL, SEX = NULL, PHENOTYPE = NULL)]

# assign APOE haplotypes
apoe[rs429358_C == 0 & rs7412_T == 0, apoe := 'apoe33']
apoe[rs429358_C == 0 & rs7412_T == 1, apoe := 'apoe23']
apoe[rs429358_C == 0 & rs7412_T == 2, apoe := 'apoe22']
apoe[rs429358_C == 1 & rs7412_T == 0, apoe := 'apoe34']
apoe[rs429358_C == 1 & rs7412_T == 1, apoe := 'apoe24'] # ambiguous w 1/3
apoe[rs429358_C == 1 & rs7412_T == 2, apoe := 'apoe12'] # rare
apoe[rs429358_C == 2 & rs7412_T == 0, apoe := 'apoe44']
apoe[rs429358_C == 2 & rs7412_T == 1, apoe := 'apoe14'] # rare
apoe[rs429358_C == 2 & rs7412_T == 2, apoe := 'apoe11'] # rare

# 


apoe <- unique(apoe)
apoe[duplicated(IID)]
apoe <- apoe[!(FID == '0' & IID %in% IID[duplicated(IID)])]
apoe <- apoe[apoe != 'apoe14']
apoe[, `:=`(rs429358_C = NULL, rs7412_T = NULL)]

apoe_covar <- merge(covar, apoe, c('FID', 'IID'))

fwrite(apoe_covar, file = "data/mega/mega_np_apoe.covar", 
       sep = ' ', quote = FALSE)
