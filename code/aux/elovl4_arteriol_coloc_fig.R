
# Part 1: Identify ELOVL4 locus
#--------------------------------
library(data.table)

results <- fread("~/nacc_rosmap.assoc.logistic")
lead_snp_pos <- results[SNP == "rs2603462", .(CHR, BP)] # lead ELOVL4 SNP
rad <- 2e5L
elovl4_locus <- results[
  CHR == lead_snp_pos$CHR[1]][
    (BP > lead_snp_pos$BP[1] - rad) & (BP < lead_snp_pos$BP[1] + rad)]

fwrite(elovl4_locus[, .(SNP)], file = "tmp/elovl4_locus.txt", 
       sep = " ", quote = FALSE, col.names = FALSE)

# Part 2: Harmonize GWAS and QTL data
#---------------------------------------
library(data.table)
library(arrow)
elovl4_rsid <- fread("tmp/elovl4_key.tmp")
elovl4_rsid[, SNP := paste0("rs", rsid)]
elovl4_id <- "ENSG00000118402.5"
qtl <- read_parquet("/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Brain_Cerebellar_Hemisphere.v8.EUR.allpairs.chr6.parquet")
setDT(qtl)
elovl4_var <- fread("tmp/elovl4_snps.tmp", header = FALSE)
qtl <- qtl[variant_id %in% elovl4_var$V1][phenotype_id == elovl4_id]
elovl4_locus <- merge(elovl4_locus, elovl4_rsid, "SNP")
elovl4_locus <- elovl4_locus[variant_id %in% qtl$variant_id]

## harmonize positions and columns for QTL and GWAS data
setorder(qtl, variant_id)
setorder(elovl4_locus, variant_id)
identical(elovl4_locus$variant_id, qtl$variant_id) # check identicalness of ids
elovl4_locus[, Phenotype := "Brain arteriolosclerosis"]
qtl[, Phenotype := "ELOVL4 expression"]
qtl[, BP := elovl4_locus$BP]
setnames(qtl, "pval_nominal", "P")
dt <- rbind(elovl4_locus[, .(Phenotype, variant_id, BP, P)],
            qtl[, .(Phenotype, variant_id, BP, P)])

# Part 3: Create colocalization plot
#-------------------------------------
library(ggplot2)

p1 <- dt[, ggplot(.SD, aes(BP / 1e6, -log10(P), color = Phenotype))]
p1 + 
  geom_point() +
  labs(x = "Chromosome 6 Position (MB, Hg19)") +
  theme_minimal()
