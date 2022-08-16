
#=====================================================================
# Intake a SAIGE null model and output residuals as phenotype in new
# phenotype-covariate data table
#=====================================================================

library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-m", "--model", help = "SAIGE null model .rda file")
parser$add_argument("--colname", help = "output column name for residuals")
parser$add_argument("--sampleID_col", default = "IID", 
                    help = "sample ID column name in covariate file")
parser$add_argument("-c", "--covar", help = "covariate file")
parser$add_argument("-o", "--output", help = "output file")
args <- parser$parse_args()

# #----------------------------------------------------
# # test args
# args$colname <- "PC1_res"
# args$covar <- "data/mega/mega_np.covar"
# args$model <- "output/multivar/vcid_PC1_stage1.rda"
# args$output <- "tmp/vcid_PC1_res.txt"
# #----------------------------------------------------

env1 <- new.env()
m <- load(args$model, env1)
nullmod <- get(m, env1)
rm(env1, m)

pheno <- data.table(v1 = nullmod$sampleID, v2 = drop(nullmod$residuals))
names(pheno) <- c(args$sampleID_col, args$colname)

covar <- fread(args$covar)
out <- merge(pheno, covar, args$sampleID_col)

fwrite(out, file = args$output, sep = " ", quote = FALSE, na = "NA")
