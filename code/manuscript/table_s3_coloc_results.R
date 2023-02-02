
library(data.table)
library(flextable)
library(magrittr)
library(stringi)
coloc <- fread("output/coloc/mega/coloc_results.txt")
coloc <- coloc[!duplicated(coloc)]

coloc_ts3_dt <- coloc[
  PPH4 >= 0.8,
  .(GWAS_Phenotype, Chromosome, QTL_Target, Tissue, QTL_Type, QTL_Source,
    P_GWAS, P_QTL, PPH4)]
coloc_ts3_dt[, Chromosome := 
               as.integer(stri_replace_first_fixed(Chromosome, "chr", ""))]
setorder(coloc_ts3_dt, GWAS_Phenotype, Chromosome, QTL_Target, Tissue)
coloc_ts3_dt[, Tissue := stri_replace_all_fixed(Tissue, "_", " ")]
coloc_ts3_dt[, P_GWAS := formatC(P_GWAS, digits = 1, format = "e")]
coloc_ts3_dt[, P_QTL := formatC(P_QTL, digits = 1, format = "e")]
coloc_ts3_dt[, PPH4 := paste0(round(PPH4 * 100, 1), "%")]
coloc_ts3_dt[PPH4 == "100%", PPH4 := ">99.9%"]
coloc_ts3_dt[, GWAS_Phenotype := 
               factor(
                 GWAS_Phenotype, 
                 labels = c("Arteriolosclerosis", "Atherosclerosis",
                            "Braak NFT Stage", "CAA", "CERAD Score", 
                            "Diffuse Amyloid Plaques", "Gross Infarction",
                            "HS", "LATE-NC", "Lewy Body", "Microinfarct"))]
setnames(coloc_ts3_dt, colnames(coloc_ts3_dt), 
         c("NPE", "Chromosome", "QTL Target", "Tissue", "QTL Type", 
           "QTL Source", "P NPE", "P QTL", "PrC"))
ts3 <- flextable(coloc_ts3_dt) %>% 
  autofit() %>% 
  padding(padding = 2, part = "all") %>%
  add_header_lines(paste0("Table S3: NPE-QTL Colocalizing Loci")) %>% 
  footnote(i = nrow_part(., part = "header"), j = 2, 
           value = 
             as_paragraph(paste0("Closest protein-coding gene according to ",
                                 "GENCODE release 40.")), 
           ref_symbols = c("a"),
           part = "header")
 ts3
 