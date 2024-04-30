library(data.table)
library(readxl)
library(stringi)
library(LDlinkR)
library(flextable)
library(magrittr)
token <- Sys.getenv("LDLINK_TOKEN")
MIN_MAF <- 0.005
MIN_R2 <- 0.97
ad_gwas_results <- 
  fread("/data_global/summary_statistics/Bellenguez/GCST90027158_buildGRCh38.tsv.gz")
ad_snps_old <- setDT(read_xlsx('raw_data/bellenguez_et_al_2022_snps.xlsx'))
ad_snps_new <- setDT(read_xlsx('raw_data/bellenguez_et_al_2022_snps.xlsx', 
                               sheet = 2))
ad_snps_old <- ad_snps_old[MAF >= MIN_MAF] # 7 SNPs MAF < 5%
ad_snps_new <- ad_snps_new[MAF >= MIN_MAF] # 5 SNPs MAF < 5%


source("code/functions/find_closest_gene.R")
source("code/functions/make_or_95_ci.R")

bim <- fread("data/adc/adc_np.bim")
setnames(bim, 
         c("V1", "V2", "V4", "V5", "V6"), 
         c("CHR", "SNP", "BP", "A1", "A2"))
# pheno_ord <- fread("data/mega/mega_np_ord.pheno")
# pheno_bin <- fread("data/mega/mega_np.pheno")
pheno = fread("shared/nacc_rosmap_act_np.pheno")

# results_files_ord <- paste0("output/gwas/mega/polmm/", 
#                             colnames(pheno_ord)[3:ncol(pheno_ord)],
#                             "_polmm_results.txt")
# results_files_bin <- paste0("output/gwas/mega/saige/",
#                             colnames(pheno_bin)[3:ncol(pheno_bin)],
#                             "_saige_results.txt")
results_files = paste0(
  "output/gwas/metal/results/", 
  colnames(pheno)[3:ncol(pheno)],
  "1.csv"
)

get_results_list <- function(files) {
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  return(file_list)
}

# collate binary and ordinal outcome results
results <- get_results_list(results_files)
for (i in 1:(ncol(pheno) - 2)) {
  results[[i]][, Phenotype := colnames(pheno)[i + 2]]
}
results <- rbindlist(results)
# results = results[, c("Chromosome", "BP", "Variant", "Allele1", "Allele2", "MAF", )]
setnames(results,
         c("chr", "bp", "MarkerName", "Freq1", "Effect", "P.value", "Allele1", "Allele2"),
         c("Chromosome", "BP", "Variant", "MAF", "Beta", "p_value", "A1", "A2"))


# format and merge with AD SNPs
ad_snps_old[, Gene := NULL]
setnames(ad_snps_old, "Known locus", "Gene")
ad_snps = rbind(ad_snps_old, ad_snps_new[, -c("Locus number")])
ad_snps[, CHR_BP := paste0(Chromosome, "_", Position)]
results[, CHR_BP := paste0(Chromosome, "_", BP)]

results_ad_snps <- results[CHR_BP %in% ad_snps$CHR_BP]

# 19 variants with MAF > 5% were not present in our GWAS, so we will query
# LDlink LDProxy to see if any reasonable proxies (LD r2 >= 97%) are available

# ad_snps_missing <- c(ad_snps_old$Variant, ad_snps_new$Variant)
# ad_snps_missing <- ad_snps[
#   !(CHR_BP %in% results_ad_snps$CHR_BP),
#   Variant
# ]
# 
# query_filename <- "combined_query_snp_list_grch38.txt"
# new_query_filename <- paste0("output/misc/", query_filename)
# 
# if (!file.exists(new_query_filename)) {
#   my_query <- LDproxy_batch(ad_snps_missing,
#                             genome_build = "grch38",
#                             append = TRUE,
#                             token = token)
#   
#   
#   file.copy(query_filename, new_query_filename, overwrite = TRUE)
#   file.remove(query_filename)
# }
# 
# 
# ldproxy <- fread(new_query_filename)
# ldproxy <- ldproxy[query_snp != RS_Number]
# ldproxy <- ldproxy[R2 >= MIN_R2 & MAF >= MIN_MAF]
# possible_ldproxies <- bim[SNP %in% ldproxy$RS_Number]
# ldproxy <- ldproxy[RS_Number %in% possible_ldproxies$SNP]
# best_proxies_index <- ldproxy[, .I[which.max(R2)], query_snp][["V1"]]
# ldproxy <- ldproxy[best_proxies_index]
# # proxies for 6 variants are available
# 
# results_ad_snps_proxy <- results[Variant %in% ldproxy$RS_Number]
# results_ad_snps_proxy <- merge(results_ad_snps_proxy,
#                                ldproxy[, .(RS_Number, query_snp, R2)],
#                                by.x = "Variant",
#                                by.y = "RS_Number")
# 
# results_ad_snps[, query_snp := Variant]
# results_ad_snps[, R2 := 1]
# results_ad_snps <- rbind(results_ad_snps, results_ad_snps_proxy)

# Use Benjamini-Hochberg false discovery rate adjustment for p-values by each
# neuropathological outcome
results_ad_snps[, Q_val := p.adjust(p_value, method = "fdr"), Phenotype]

# merge with Bellenguez et al 2022 stage 1 summary statistics and check for
# effect direction concordance
ad_snps_used = ad_snps[(CHR_BP %in% results_ad_snps$CHR_BP)]
ad_gwas_results <- ad_gwas_results[variant_id %in% ad_snps_used$Variant]
ad_gwas_results[, CHR_BP := paste0(chromosome, "_", base_pair_location)]
results_ad_snps <- merge(
  results_ad_snps,
  ad_gwas_results[, .(variant_id, effect_allele, other_allele, beta, p_value, CHR_BP)],
  by= "CHR_BP"
)

# For both studies, does each SNPs either have same effect/other allele 
# or switched alleles? Yes

# number of variants where alleles don't match
results_ad_snps[A1 != effect_allele | A2 != other_allele, .N] # 250

# number of variants where alleles match after switching minor and major alleles
results_ad_snps[
  A1 != effect_allele | A2 != other_allele
][
  A2 == effect_allele & A1 == other_allele, .N] # 240

# remove variants that still don't match after swiching alleles
remove_variants = results_ad_snps[
  A1 != effect_allele | A2 != other_allele
][
  A2 != effect_allele | A1 != other_allele,
  unique(Variant)
]

results_ad_snps = results_ad_snps[!(Variant %in% remove_variants)]

# reverse effect sign on variants where alleles were switched
results_ad_snps[A2 == effect_allele & A1 == other_allele,
                beta := -beta]

# create concordance indicator
results_ad_snps[, Concordant := beta * Beta > 0]
results_ad_snps[, Concordant := factor(Concordant, labels = c("No", "Yes"))]

# format columns for output table
results_ad_snps[, `Effect/other allele` := paste0(effect_allele, "/", other_allele)]
results_ad_snps[, OR := round(exp(Beta), 2)]
results_ad_snps[, or := round(exp(beta), 2)]
results_ad_snps[, beta := round(beta, 2)]
results_ad_snps[, Beta := round(Beta, 2)]
results_ad_snps[, P := formatC(p_value.x, digits = 2)]
results_ad_snps[, p_value := formatC(as.numeric(p_value.y), digits = 2)]
results_ad_snps[, Q_val := formatC(Q_val, digits = 2)]
results_ad_snps[, MAF := round(MAF, 2)]
# add gene info from Bellenguez tables
results_ad_snps = merge(results_ad_snps, ad_snps[, .(CHR_BP, Gene)], "CHR_BP")

results_ad_snps <- results_ad_snps[
  ,
  .(Phenotype, Chromosome, Gene, BP, variant_id, MAF, 
    `Effect/other allele`, Beta, OR, P, Q_val, beta, or, p_value, Concordant)
]
setorder(results_ad_snps, Phenotype, Chromosome, BP)
results_ad_snps[
  , 
  Phenotype := factor(
    Phenotype, 
    labels = c("Arteriolosclerosis", "Atherosclerosis", "Braak NFT Stage", 
               "CAA", "CERAD Score", "Amyloid-Beta Plaques",
               "Gross Infarction", "HS", "LATE-NC", "Lewy Body", 
               "Microinfarct"
    )
  )
]

setnames(results_ad_snps, 
         colnames(results_ad_snps),
         c("NPE", "Chromosome", "Locus", "Position", "Variant", 
           "EAF", "Effect/Other Allele", 
           "NPE Beta", "NPE OR", "NPE P-value", "NPE Q-value", "ADD Beta", 
           "ADD OR", "ADD P-value", "NPE-ADD Concordant Effect Direction"
           )
)


fwrite(results_ad_snps,
       file = "doc/Supplementary_Table_Bellenguez_AD_GWAS_SNPs.csv",
       quote = FALSE)

results_ad_snps_signif <- results_ad_snps[as.numeric(`NPE Q-value`) < 0.05]

fwrite(results_ad_snps_signif,
       file = "doc/Supplementary_Table_Bellenguez_AD_GWAS_SNPs_FDR0.05.csv",
       quote = FALSE)

# format table of variants with Q values < 0.05 for main text
results_ad_snps_table <- results_ad_snps_signif %>% 
  flextable() %>% 
  autofit() %>% 
  padding(padding = 2, part = "all") %>%
  add_header_lines(paste0("Table 2: Associations between NPE and known ADD loci"))%>% 
  style(j = 3, 
        pr_t = fp_text_default(
          italic = TRUE)) %>% 
  add_footer_lines(
    value = as_paragraph(
      paste0(
        "ADD = Alzheimer's Disease and related dementias; ",
        "EAF = Effect allele frequency. ",
        "NPE P-values, betas, and OR are from meta-analysis. ",
        "NPE Q-values are produced by applying Benjamini-Hochberg adjustments ",
        "for each endophenotype separately. ",
        "ADD P-values, betas, and OR are from Bellenguez et al (2022) stage I GWAS (N = 487,511). ",
        "OR are with respect to the Bellenguez effect allele. "
      )
    )
  ) %>%
  footnote(i = nrow_part(., part = "header"), j = 3,
           value =
             as_paragraph(
               paste0(
                 "Either known locus or closest protein-coding gene according ",
                      "to GENCODE release 40.")),
           ref_symbols = c("a"),
           part = "header",
           inline = TRUE,
           sep = " ") %>%
  footnote(i = nrow_part(., part = "header"), j = 4,
           value =
             as_paragraph(
               paste0(
                 "Position of lead variant using GRCh38 assembly.")),
           ref_symbols = c("b"),
           part = "header",
           inline = TRUE,
           sep = " ") %>% 
  footnote(i = nrow_part(., part = "header"), j = 6,
           value =
             as_paragraph(
               paste0(
                 "Effect allele frequency in NPE meta-analysis.")),
           ref_symbols = c("c"),
           part = "header",
           inline = TRUE,
           sep = " ")
  

save(results_ad_snps_table, 
     file = "output/manuscript/results_ad_snps_table.Rdata")
