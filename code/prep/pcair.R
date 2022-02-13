#! /usr/bin/Rscript --vanilla

#==============================================================
# perform PC-AiR 
#==============================================================

library(data.table)
library(GENESIS)
library(SNPRelate)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-b", "--bfile",
                    help="PLINK binary file set prefix")
parser$add_argument("-k", "--kinship",
                    help="KING kinship matrix .kin file path")
parser$add_argument("-o", "--out",
                    help="Output .Rdata file path")
parser$add_argument("-g", "--gdsfile", default="tmp/pcair.gds",
                    help="path to create GDS file")
args <- parser$parse_args()
#-----------
# PC-AiR
#-----------
# mostly taken from here: 
# https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html
gdsfile <- args$gdsfile
plink_path <- args$bfile
snps <- fread(paste0(plink_path, ".bim"))
snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"),
              bim.fn = paste0(plink_path, ".bim"),
              fam.fn = paste0(plink_path, ".fam"),
              out.gdsfn = gdsfile)

# create kinship matrix
KINGmat <- kingToMatrix(args$kinship, estimator = "Kinship")

geno <- GWASTools::GdsGenotypeReader(filename = gdsfile)

mypcair <- pcair(geno, kinobj = KINGmat, divobj = KINGmat, snp.include = snps$V2)
save(mypcair, file = args$out)
