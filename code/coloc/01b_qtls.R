##-----------------------------------
## Tidy QTL files for input SNPs
##-----------------------------------
library(pacman)
p_load(data.table, magrittr, stringi)
qtls <- fread(commandArgs(trailingOnly = T)[1], header = F)
setnames(qtls, colnames(qtls), 
         c("phenotype_id", "variant_id",	"tss_distance",	"maf",	"ma_samples",
           "ma_count",	"pval_nominal",	"slope",	"slope_se",	"pval_nominal_threshold",
           "min_pval_nominal", "pval_beta"))

# split the file:phenotype string into separate columns
genes_files <- matrix(unlist(stri_split_fixed(qtls$phenotype_id, ":", n=2)), ncol=2, byrow = T)
qtls[, filename := genes_files[, 1]]
qtls[, phenotype_id := genes_files[, 2]]

# create variable for QTL type
qtls[grep("eQTL", filename), qtl_type := "eQTL"]
qtls[grep("sQTL", filename), qtl_type := "sQTL"]

# create a variable for tissue type (eQTLs and sQTLs have different file prefixes)
qtls[qtl_type == "eQTL", tissue := stri_replace_first_regex(
  filename, "/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_Eur/eqtls/", "")]
qtls[qtl_type == "sQTL", tissue := stri_replace_first_regex(
  filename, "/data_global/GTEx/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_sQTL_Eur/GTEx_Analysis_v8_sQTL_EUR/", "")]
qtls[, tissue := stri_replace_last_fixed(tissue, ".v8.EUR.signif_pairs.txt", "")]

# create gene variable
qtls[qtl_type == "eQTL", gene := phenotype_id]
qtls[qtl_type == "sQTL", gene := stri_replace_all_regex(phenotype_id, ".*ENSG", "ENSG")]

# merge with rsid_key
rsid_key <- fread(commandArgs(trailingOnly = T)[2])
qtls <- merge(qtls, rsid_key, "variant_id")
setcolorder(qtls, c("rsid", "variant_id", "phenotype_id", "gene", "tissue", "pval_nominal", "slope", "slope_se"))
setorder(qtls, rsid, tissue)

setkey(qtls, phenotype_id, tissue, pval_nominal)

top_qtl_pvals <- qtls[, min(pval_nominal), .(phenotype_id, tissue)] %>%
  setnames("V1", "pval_nominal")

top_qtls <- merge(qtls, top_qtl_pvals, colnames(top_qtl_pvals)) %>% 
  .[!duplicated(.[, .(phenotype_id, tissue, pval_nominal)])]

# create list of chromosomes needed for colocalization
top_qtls[, chr := 
           stri_replace_first_regex(variant_id, "_.*", "") %>% 
           stri_replace_first_fixed(., "chr", "") %>% 
           as.integer() 
           
]

setcolorder(top_qtls, c("phenotype_id", "tissue", "chr", "qtl_type"))

chrs <- top_qtls[, chr] %>% 
  unique() %>% 
  as.data.table()

# write files
fwrite(top_qtls, file = "data/tmp/top_qtls.tmp", sep = " ", row.names = F, col.names = T, quote = F)
fwrite(chrs, file = "data/tmp/chrs.tmp", sep=" ", row.names = F, col.names = F, quote = F)
