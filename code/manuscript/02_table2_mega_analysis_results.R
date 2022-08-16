#==========================================
# Produce Table 2: GWAS results summary
#==========================================

library(data.table)
library(magrittr)
library(rtracklayer)
library(flextable)
source("code/functions/find_closest_gene.R")
source("code/functions/make_or_95_ci.R")

bim <- fread("data/mega/mega_np.bim")
setnames(bim, 
         c("V1", "V2", "V4", "V5", "V6"), 
         c("CHR", "SNP", "BP", "A1", "A2"))
pheno_ord <- fread("data/mega/mega_np_ord.pheno")
pheno_bin <- fread("data/mega/mega_np.pheno")

results_files_ord <- paste0("output/gwas/mega/polmm/", 
                            colnames(pheno_ord)[3:ncol(pheno_ord)],
                            "_polmm_results.txt")
results_files_bin <- paste0("output/gwas/mega/saige/",
                            colnames(pheno_bin)[3:ncol(pheno_bin)],
                            "_saige_results.txt")
clump_files_ord <- paste0("output/gwas/mega/polmm/", 
                          colnames(pheno_ord)[3:ncol(pheno_ord)],
                          "_polmm_results.clumped")
clump_files_bin <- paste0("output/gwas/mega/saige/",
                          colnames(pheno_bin)[3:ncol(pheno_bin)],
                          "_saige_results.clumped")

get_results_list <- function(files) {
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  return(file_list)
}

# collate binary and ordinal outcome results
results_ord <- get_results_list(results_files_ord)
for (i in 1:(ncol(pheno_ord) - 2)) {
  results_ord[[i]][, Phenotype := colnames(pheno_ord)[i + 2]]
}
results_ord <- rbindlist(results_ord)

results_bin <- get_results_list(results_files_bin)
for (i in 1:(ncol(pheno_bin) - 2)) {
  results_bin[[i]][, Phenotype := colnames(pheno_bin)[i + 2]]
}
results_bin <- rbindlist(results_bin)

# add significant locus from APOE-adjusted CAA analysis
caa_apoe_results <- 
  fread("output/gwas/mega/polmm/caa_ord_apoe_polmm_results.txt")

# collate clumped results
results_clump_ord <- get_results_list(clump_files_ord)
names(results_clump_ord) <- colnames(pheno_ord)[3:ncol(pheno_ord)]
for (i in 1:(ncol(pheno_ord) - 2)) {
  results_clump_ord[[i]][, Phenotype := colnames(pheno_ord)[i + 2]]
}
results_clump_ord <- rbindlist(results_clump_ord)
results_ord_top <- merge(results_ord, 
                         results_clump_ord[, .(SNP, Phenotype)], 
                         by.x = c("SNPID", "Phenotype"),
                         by.y = c("SNP", "Phenotype"))
results_ord_top <- 
  results_ord_top[, .(Phenotype, SNPID, chr, MAF, beta, pval.spa)]
setnames(results_ord_top, 
         colnames(results_ord_top),
         c("Phenotype", "Variant", "Chromosome", "MAF", "Beta", "P"))

results_clump_bin <- get_results_list(clump_files_bin)
for (i in 1:(ncol(pheno_bin) - 2)) {
  results_clump_bin[[i]][, Phenotype := colnames(pheno_bin)[i + 2]]
}
results_clump_bin <- rbindlist(results_clump_bin)
results_bin_top <- merge(results_bin, 
                         results_clump_bin[, .(SNP, Phenotype)], 
                         by.x = c("MarkerID", "Phenotype"),
                         by.y = c("SNP", "Phenotype"))

results_bin_top <- 
  results_bin_top[, .(Phenotype, MarkerID, CHR, 1-AF_Allele2, 
                      BETA, p.value)]
setnames(results_bin_top, 
         colnames(results_bin_top),
         c("Phenotype", "Variant", "Chromosome", "MAF", "Beta", "P"))

caa_apoe_clumped <- 
  fread("output/gwas/mega/polmm/caa_ord_apoe_polmm_results.clumped")
caa_apoe_top <- merge(caa_apoe_results, 
                      caa_apoe_clumped[, .(SNP)], 
                      by.x = c("SNPID"),
                      by.y = c("SNP"))
caa_apoe_top[, Phenotype := "caa_ord"]
caa_apoe_top <- 
  caa_apoe_top[, .(Phenotype, SNPID, chr, MAF, beta, pval.spa)]
setnames(caa_apoe_top, 
         colnames(caa_apoe_top),
         c("Phenotype", "Variant", "Chromosome", "MAF", "Beta", "P"))


# combine all top results
results_top <- rbindlist(list(results_ord_top, 
                              results_bin_top, 
                              caa_apoe_top))

# merge top results with .bim file
results_top <- merge(results_top, 
                     bim, 
                     by.x = c("Chromosome", "Variant"), 
                     by.y = c("CHR", "SNP"))

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
results_top[, Phenotype := 
              factor(
                Phenotype, 
                labels = c("Arteriolosclerosis", "Atherosclerosis",
                           "Braak Stage", "CAA", "CERAD Score", 
                           "Diffuse Amyloid Plaques", "Gross Infarction",
                           "HS", "LATE-NC", "Lewy Body", "Microinfarct",
                           "PART", "VCID"
                ))]
results_top[, `Minor/major allele` := paste0(A1, "/", A2)]
results_top[, `95% CI` := make_95CI(OR, P)]
setnames(results_top, "BP", "Position")
setcolorder(results_top, c("Chromosome", "Position", "Phenotype"))
results_signif <- 
  results_top[P < 5e-8, 
              .(Phenotype, Variant, Chromosome, Position, Gene, 
                `Minor/major allele`, OR, `95% CI`, P)
  ]
results_signif[, P := formatC(P, digits = 2)]

# for phenotypes with more than one hit in APOE region, 
# remove all but top hits
results_signif <- results_signif[
  !((Phenotype == "CAA" & Variant %in% c("rs111997200", "rs283810")) |
      (Phenotype == "Braak Stage" & Variant == "rs11668327") |
      (Phenotype == "CERAD Score" & Variant == "rs365653"))
]
setorder(results_signif, Chromosome, Position, Phenotype)

t2 <- flextable(results_signif) %>% 
  autofit() %>% 
  padding(padding = 2, part = "all") %>%
  add_header_lines(paste0("Table 2: Significant NPE-Associated Loci in ", 
                          "Mega-Analysis")) %>% 
  # add_footer_lines(paste0("Gene indicated closes protein-coding gene using ",
  #                         "Gencode v40. ORs are with respect to minor allele.",
  #                         " Genome positions are based on build HG38.")) %>% 
  # add_footer_lines(paste0("Key: OR, odds ratio; CI, confidence interval; ",
  #                         "LATE-NC, LATE neuropathologic change; HS, ",
  #                         "hippocampal sclerosis; CAA, cerebral amyloid ",
  #                         "angiopathy.")) %>% 
  footnote(i = nrow_part(., part = "header"), j = 4, 
           value = 
             as_paragraph("Genome positions are based on build HG38."), 
           ref_symbols = c("a"),
           part = "header") %>% 
  footnote(i = nrow_part(., part = "header"), j = 5, 
           value = 
             as_paragraph(paste0("Closest protein-coding gene according to ",
                                 "GENCODE release 40.")), 
           ref_symbols = c("b"),
           part = "header") %>% 
  footnote(i = nrow_part(., part = "header"), j = 7, 
           value = 
             as_paragraph(paste0("ORs are with respect to minor allele.")), 
           ref_symbols = c("c"),
           part = "header") %>% 
  footnote(i = nrow_part(., part = "body"), j = 1, 
           value = 
             as_paragraph("Result from APOE diplotype-adjusted analysis"), 
           ref_symbols = c("d")) %>% 
  footnote(i = (nrow_part(., part = "body")-5):nrow_part(., part = "body"), 
           j = 5, 
           value = 
             as_paragraph("Locus in APOE region"), 
           ref_symbols = c("e"))

t2
