#===============================================
# Create data set with IIDs and PC-AiR PCs 1-5
#===============================================
# mostly taken from here: 
# https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html

library(pacman)
p_load(data.table, GENESIS, SNPRelate)
cargs <- commandArgs(trailingOnly = TRUE)
plink_prefix <- cargs[1]
kinship_file <- cargs[2]

plink_prefix_pruned <- paste0(plink_prefix, "_pruned")

gdsfile <- paste0("data/tmp/", basename(plink_prefix), "_pcair.gds")
snps <- fread(paste0(plink_prefix_pruned, ".bim"))
snpgdsBED2GDS(bed.fn = paste0(plink_prefix_pruned, ".bed"),
              bim.fn = paste0(plink_prefix_pruned, ".bim"),
              fam.fn = paste0(plink_prefix_pruned, ".fam"),
              out.gdsfn = gdsfile)

# create kinship matrix
KINGmat <- kingToMatrix(kinship_file, estimator = "Kinship")

geno <- GWASTools::GdsGenotypeReader(filename = gdsfile)

mypcair <- pcair(geno, 
                 kinobj = KINGmat, 
                 divobj = KINGmat, 
                 snp.include = snps$V2)

pcs <- as.data.table(mypcair$vectors[, 1:5], keep.rownames = TRUE)
setnames(pcs, colnames(pcs), c("IID", "PC1", "PC2", "PC3", "PC4", "PC5"))

file.remove(gdsfile)

fwrite(pcs, 
       file = paste0(plink_prefix, "_pcs.txt"),
       sep = " ",
       quote = FALSE)