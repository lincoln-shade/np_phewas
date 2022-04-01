##-----------------------------------
## Tidy QTL files for input SNPs
##-----------------------------------

library(data.table)
library(magrittr)
library(stringi)
args <- commandArgs(trailingOnly = TRUE)
phenotype <- args[3]
qtls <- fread(args[1], header = F)
setnames(qtls, colnames(qtls), 
         c("filename", "phenotype_id", "variant_id",	"tss_distance",	"maf",
           "ma_samples", "ma_count",	"pval_nominal",	"slope",	"slope_se",	
           "pval_nominal_threshold", "min_pval_nominal", "pval_beta"))

# create variable for QTL type
qtls[!(grep("clu", phenotype_id)), qtl_type := "eQTL"]
qtls[grep("clu", phenotype_id), qtl_type := "sQTL"]

qtls[, tissue := stri_replace_last_fixed(filename, ".v8.EUR.signif_pairs.txt", "")]

# create gene variable
qtls[qtl_type == "eQTL", gene := phenotype_id]
qtls[qtl_type == "sQTL", gene := 
       stri_replace_all_regex(phenotype_id, ".*ENSG", "ENSG")]

# merge with rsid_key
rsid_key <- fread(args[2])
qtls <- merge(qtls, rsid_key, "variant_id")
setcolorder(qtls, c("rsid", "variant_id", "phenotype_id", "gene", "tissue", 
                    "pval_nominal", "slope", "slope_se"))
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
fwrite(top_qtls, 
       file = paste0("tmp/top_qtls_", phenotype, ".tmp"), 
       sep = " ", 
       col.names = T, 
       quote = F
)

fwrite(chrs, 
       file = paste0("tmp/chrs_", phenotype, ".tmp"), 
       sep=" ", 
       col.names = F, 
       quote = F
)
