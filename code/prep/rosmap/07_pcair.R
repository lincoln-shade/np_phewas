#==============================================================
# perform PC-AiR 
#==============================================================


pacman::p_load(data.table, GENESIS, SNPRelate)
study <- 'rosmap'
#-----------
# PC-AiR
#-----------
# mostly taken from here: 
# https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html
gdsfile <- "tmp/pcair.gds"
plink_path <- "data/rosmap/rosmap_np_pruned"
snps <- fread(paste0(plink_path, ".bim"))
snpgdsBED2GDS(bed.fn = paste0(plink_path, ".bed"),
              bim.fn = paste0(plink_path, ".bim"),
              fam.fn = paste0(plink_path, ".fam"),
              out.gdsfn = gdsfile)

# create kinship matrix
KINGmat <- kingToMatrix("data/rosmap/rosmap_np.kin", estimator = "Kinship")

geno <- GWASTools::GdsGenotypeReader(filename = gdsfile)

mypcair <- pcair(geno, kinobj = KINGmat, divobj = KINGmat, snp.include = snps$V2)
save(mypcair, file = "data/rosmap/mypcair.Rdata")
