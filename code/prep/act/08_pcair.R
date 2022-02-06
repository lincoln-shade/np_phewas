#===============================================
# Perform PCAiR
#===============================================
# mostly taken from here: 
# https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html

library(pacman)
p_load(data.table, GENESIS, SNPRelate)
cargs <- commandArgs(trailingOnly = TRUE)
plink_prefix <- 'data/act/act_np'
kinship_file <- 'data/act/act_np.kin'

plink_prefix_pruned <- paste0(plink_prefix, "_pruned")

gdsfile <- paste0("tmp/", basename(plink_prefix), "_pcair.gds")
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

save(mypcair, file = paste0('data/act/mypcair.Rdata'))

file.remove(gdsfile)