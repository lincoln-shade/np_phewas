#==========================================
# Produce Table 2: GWAS results summary
#==========================================

library(data.table)
library(magrittr)
library(rtracklayer)
library(flextable)
source("code/functions/find_closest_gene.R")
source("code/functions/make_or_95_ci.R")

bim <- fread("data/adc/adc_np.bim")
setnames(bim, 
         c("V1", "V2", "V4", "V5", "V6"), 
         c("CHR", "SNP", "BP", "A1", "A2"))

pheno = fread("shared/nacc_rosmap_act_np.pheno")

results_files <- paste0(
  "output/gwas/metal/results/", colnames(pheno)[3:ncol(pheno)], "1.csv"
)

clump_files <- paste0(
  "output/gwas/metal/results/", colnames(pheno)[3:ncol(pheno)], "1_clump.clumped"
)

get_results_list <- function(files) {
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  return(file_list)
}

# collate results
results <- get_results_list(results_files)

for (i in 1:(ncol(pheno) - 2)) {
  results[[i]][, Phenotype := colnames(pheno)[i + 2]]
}
results <- rbindlist(results)


# add significant locus from APOE-adjusted CAA analysis
caa_apoe_results <- fread("output/gwas/metal/results/caa_apoe1.csv")

# collate clumped results
results_clump <- get_results_list(clump_files)

names(results_clump) <- colnames(pheno)[3:ncol(pheno)]

for (i in 1:(ncol(pheno) - 2)) {
  results_clump[[i]][, Phenotype := colnames(pheno)[i + 2]]
}

results_clump <- rbindlist(results_clump)

results_top <- merge(
  results, 
  results_clump[, .(SNP, Phenotype)], 
  by.x = c("MarkerName", "Phenotype"),
  by.y = c("SNP", "Phenotype")
)

results_top <- 
  results_top[, .(Phenotype, MarkerName, chr, bp, Freq1, Allele1, Allele2, Effect, StdErr, P.value)]
setnames(results_top, 
         colnames(results_top),
         c("Phenotype", "Variant", "Chromosome", "BP", "MAF", "A1", "A2", "Beta", "SE", "P"))

caa_apoe_top <- caa_apoe_results[P.value == min(P.value)]
caa_apoe_top[, Phenotype := "caa"]

caa_apoe_top <- 
  caa_apoe_top[, .(Phenotype, MarkerName, chr, bp, Freq1, Allele1, Allele2, Effect, StdErr, P.value)]

setnames(
  caa_apoe_top, 
  colnames(caa_apoe_top),
  c("Phenotype", "Variant", "Chromosome", "BP", "MAF", "A1", "A2", "Beta", "SE", "P")
)


# combine all top results
results_top <- rbindlist(list(results_top, caa_apoe_top))

# # merge top results with .bim file
# results_top <- merge(results_top, 
#                      bim, 
#                      by.x = c("Chromosome", "Variant"), 
#                      by.y = c("CHR", "SNP"))

# ensure minor allele is effect allele
results_top[MAF > 0.5, Beta := -Beta]
majora = results_top[MAF > 0.5, A1]
minora = results_top[MAF > 0.5, A2]
results_top[MAF > 0.5, A1 := A2]
results_top[MAF > 0.5, A2 := ..majora]
results_top[MAF > 0.5, MAF := 1-MAF]

# import list of gene names and positions and 
# add closest protein-coding gene to table
genes <- 
  rtracklayer::import("~/downloads/gencode.v40.annotation.gtf.gz",
                      feature.type = "gene") %>% 
  as.data.table() %>% 
  .[gene_type == "protein_coding"]
genes <- genes[, .(gene_name, seqnames, start, end)] %>% 
  melt(., id.vars = c("gene_name", "seqnames"), 
       measure.vars = c("start", "end"))
genes[, chr := as.integer(seqnames)]
results_top[, Gene := find_closest_gene(Chromosome, BP, genes)]

results_top[, OR := round(exp(Beta), 2)]
results_top[, P := signif(P, 2)]

results_top[
  , 
  Phenotype := factor(
    Phenotype, 
    labels = c("Arteriolosclerosis", "Atherosclerosis",
               "Braak NFT Stage", "CAA", "CERAD Score", "Gross Infarction",
               "HS", "LATE-NC", "Lewy Body", "Microinfarct")
  )
]

results_top[, `Minor/major allele` := paste0(A1, "/", A2)]
results_top[, `95% CI` := make_or_95_ci(OR, round(exp(Beta - 1.96 * SE), 2), round(exp(Beta + 1.96 * SE), 2), P)]
setnames(results_top, "BP", "Position")
setcolorder(results_top, c("Chromosome", "Position", "Phenotype"))
results_signif <- 
  results_top[P < 1e-7, 
              .(Phenotype, Gene, Variant, Chromosome, Position, 
                `Minor/major allele`, OR, `95% CI`, P)
  ]
results_signif[, P := formatC(P, digits = 2)]

# for phenotypes with more than one hit in APOE region, 
# remove all but top hits
results_signif <- results_signif[
  !(
    (Phenotype == "CAA" & Variant %in% c("chr19:44860563:T:G", "rs283810", "rs3760625")) |
      (Phenotype == "Braak NFT Stage" & Variant %in% c("chr19:44909976:G:T", "rs3852860")) |
      (Phenotype == "CERAD Score" & Variant %in% c("chr19:44909976:G:T", "rs3852860")))
]
setorder(results_signif, Chromosome, Position, Phenotype)

t2 <- flextable(results_signif) %>% 
  autofit() %>% 
  padding(padding = 2, part = "all") %>%
  add_header_lines(paste0("Table 1: Significant NPE-Associated Loci in ", 
                          "Stage 3 pooled GWAS")) %>% 
  footnote(i = nrow_part(., part = "header"), j = 2, 
           value = 
             as_paragraph(paste0("Closest protein-coding gene according to ",
                                 "GENCODE release 40. ")), 
           ref_symbols = c("a"),
           part = "header", 
           inline = TRUE,
           sep = " ") %>% 
  footnote(i = nrow_part(., part = "header"), j = 5, 
           value = 
             as_paragraph("Genome positions are based on build HG38."), 
           ref_symbols = c("b"),
           part = "header", 
           inline = TRUE,
           sep = " ") %>% 
  footnote(i = nrow_part(., part = "header"), j = 7, 
           value = 
             as_paragraph(paste0("ORs are with respect to minor allele.")), 
           ref_symbols = c("c"),
           part = "header", 
           inline = TRUE,
           sep = " ") %>% 
  footnote(i = nrow_part(., part = "body"), j = 1, 
           value = 
             as_paragraph("Result from APOE diplotype-adjusted analysis."), 
           ref_symbols = c("d"),
           inline = TRUE,
           sep = " ") %>% 
  # footnote(i = (nrow_part(., part = "body")-5):nrow_part(., part = "body"), 
  #          j = 2, 
  #          value = 
  #            as_paragraph("Locus in APOE region"), 
  #          ref_symbols = c("e"), 
  #          inline = TRUE,
  #          sep = " ") %>% 
  style(j = 2, 
       pr_t = fp_text_default(
         italic = TRUE)) 

save(t2, file = "output/manuscript/table2_meta_analysis_results.Rdata")
