library(GRAB)
library(data.table)
library(stringi)
source("code/functions/strip_alleles.R")
PhenoData = fread("data/act/act_np.pheno", na.strings = "-1")
covar_data = fread("data/act/act_np.covar")
PhenoData = merge(PhenoData, covar_data, c("FID", "IID"))
PhenoData[, diffuse_abeta := factor(diffuse_abeta, levels = c(0, 1, 2, 3))]
GenoFile = "data/act/act_np_pruned.bed"
getTempFilesFullGRM("data/act/act_np_pruned", 
                      nPartsGRM = 1, 
                      partParallel = 1)

getSparseGRM("data/act/act_np_pruned", 
             nPartsGRM = 1,
             SparseGRMFile = "tmp/sparse_grm.txt",
             relatednessCutoff = 0.01)

obj.POLMM = GRAB.NullModel(formula = diffuse_abeta ~ msex + age_death + ACT2 + ACT3 + pc1 + pc2 +pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                           data = PhenoData, 
                           subjData = PhenoData$IID, 
                           method = "POLMM", 
                           traitType = "ordinal",
                           GenoFile = GenoFile,
                           SparseGRMFile =  "tmp/sparse_grm.txt",
                           control = list(showInfo = FALSE, 
                                          LOCO = FALSE))

gwas_plink = "data/act/act_np.bed"
OutputFile = "output/gwas/act/polmm/diffuse_abeta_polmm_results.txt"
GRAB.Marker(
  obj.POLMM, 
  GenoFile = gwas_plink,
  OutputFile = OutputFile
)
results = fread("output/gwas/act/polmm/diffuse_abeta_polmm_results.txt")
snp_info = stri_split_fixed(results$Info, ":", simplify = T)
results[, chr := snp_info[, 1]]
results[, A2 := snp_info[, 3]]
results[, A1 := snp_info[, 4]]
results[, N := obj.POLMM$N]
setnames(results, c("Marker", "beta", "seBeta", "Pvalue", "AltFreq"), c("SNPID", "beta", "SE", "pval.spa", "MAF"))
results = results[!is.na(pval.spa)]
results[grep("rs", SNPID), SNPID := strip_alleles(SNPID)]
fwrite(results, file = "output/gwas/metal/input/act_diffuse_abeta.csv")
