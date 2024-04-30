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
pheno_ord <- fread("data/mega/nacc_np_ord.pheno")
pheno_bin <- fread("data/mega/nacc_np.pheno")

#-------------------
# nacc results
#--------------------
results_files_ord <- paste0("output/gwas/adc/polmm/", 
                            colnames(pheno_ord)[3:ncol(pheno_ord)],
                            "_polmm_results.txt")
results_files_bin <- paste0("output/gwas/adc/saige/",
                            colnames(pheno_bin)[3:ncol(pheno_bin)],
                            "_saige_results.txt")
clump_files_ord <- paste0("output/gwas/adc/polmm/", 
                          colnames(pheno_ord)[3:ncol(pheno_ord)],
                          "_polmm_results.clumped")
clump_files_bin <- paste0("output/gwas/adc/saige/",
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
# In SAIGE, effect sizes are wrt major allele (allele 2), so need to
# reverse sign of beta to get effect sizes wrt to minor allele
results_bin_top[MAF < 0.5, Beta := -Beta]


# combine all top results
results_top <- rbindlist(list(results_ord_top, 
                              results_bin_top))

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

results_top[, `Minor/major allele` := paste0(A1, "/", A2)]
results_top[, `95% CI` := make_95CI(OR, P)]
setnames(results_top, "BP", "Position")
setcolorder(results_top, c("Chromosome", "Position", "Phenotype"))
results_signif <- 
  results_top[P < 1e-5, 
              .(Phenotype, Gene, Variant, Chromosome, Position, 
                `Minor/major allele`, MAF, OR, `95% CI`, P)
  ]
results_signif[, P := formatC(P, digits = 2)]
results_signif[, MAF := round(MAF, 3)]
setorder(results_signif, Phenotype, Chromosome, Position)

#--------------------------------------------------
# ROSMAP results
#--------------------------------------------------
results_files_ord <- paste0("output/gwas/rosmap/polmm/", 
                            colnames(pheno_ord)[3:ncol(pheno_ord)],
                            "_polmm_results.txt")
results_files_bin <- paste0("output/gwas/rosmap/saige/",
                            colnames(pheno_bin)[3:ncol(pheno_bin)],
                            "_saige_results.txt")

# collate binary and ordinal outcome results
results_ord <- get_results_list(results_files_ord)
for (i in 1:(ncol(pheno_ord) - 2)) {
  results_ord[[i]][, Phenotype := colnames(pheno_ord)[i + 2]]
}
results_ord <- rbindlist(results_ord)
results_ord[, OR_rosmap := round(exp(beta), 2)]
results_ord[, P_rosmap := signif(pval.spa, 2)]
setnames(results_ord, "SNPID", "Variant")

results_bin <- get_results_list(results_files_bin)
for (i in 1:(ncol(pheno_bin) - 2)) {
  phenotype <- colnames(pheno_bin)[i + 2]
  results_bin[[i]][, Phenotype := phenotype]
  results_bin[[i]] <- results_bin[[i]][
    MarkerID %in% results_signif[Phenotype == phenotype, Variant]]
}
results_bin <- rbindlist(results_bin)
setnames(results_bin, c("CHR", "MarkerID"), c("Chromosome", "Variant"))
results_bin[, BETA := -BETA]
results_bin[, OR_rosmap := round(exp(BETA), 2)]
results_bin[, P_rosmap := signif(p.value, 2)]

rbind_cols <- c("Phenotype", "Variant", "OR_rosmap", "P_rosmap")
results_rosmap <- rbind(results_bin[, ..rbind_cols],
                        results_ord[, ..rbind_cols])

#-------------------------------------------------
# ACT results
#-------------------------------------------------
pheno_ord_no_late <- pheno_ord[, -c("late")] # no LATE data in ACT... yet
results_files_ord <- paste0("output/gwas/act/polmm/", 
                            colnames(pheno_ord_no_late)[3:ncol(pheno_ord_no_late)],
                            "_polmm_results.txt")
results_files_bin <- paste0("output/gwas/act/saige/",
                            colnames(pheno_bin)[3:ncol(pheno_bin)],
                            "_saige_results.txt")

# collate binary and ordinal outcome results
results_ord <- get_results_list(results_files_ord)

for (i in 1:(ncol(pheno_ord_no_late) - 2)) {
  results_ord[[i]][, Phenotype := colnames(pheno_ord_no_late)[i + 2]]
}
results_ord <- rbindlist(results_ord)
results_ord[, OR_act := round(exp(beta), 2)]
results_ord[, P_act := signif(pval.spa, 2)]
setnames(results_ord, "SNPID", "Variant")

results_bin <- get_results_list(results_files_bin)
for (i in 1:(ncol(pheno_bin) - 2)) {
  phenotype <- colnames(pheno_bin)[i + 2]
  results_bin[[i]][, Phenotype := phenotype]
  results_bin[[i]] <- results_bin[[i]][
    MarkerID %in% results_signif[Phenotype == phenotype, Variant]]
}
results_bin <- rbindlist(results_bin)
setnames(results_bin, c("CHR", "MarkerID"), c("Chromosome", "Variant"))
results_bin[, BETA := -BETA]
results_bin[, OR_act := round(exp(BETA), 2)]
results_bin[, P_act := signif(p.value, 2)]

rbind_cols <- c("Phenotype", "Variant", "OR_act", "P_act")
results_act <- rbind(results_bin[, ..rbind_cols],
                        results_ord[, ..rbind_cols])

#-------------------------------------------------
# ADNI results
#-------------------------------------------------
pheno_ord_adni <- pheno_ord[, -c("diffuse_abeta")]
results_files_ord <- paste0("output/gwas/adni/polmm/", 
                            colnames(pheno_ord_adni)[3:ncol(pheno_ord_adni)],
                            "_polmm_results.txt")
results_files_bin <- paste0("output/gwas/adni/saige/",
                            colnames(pheno_bin)[3:ncol(pheno_bin)],
                            "_saige_results.txt")

# collate binary and ordinal outcome results
results_ord <- get_results_list(results_files_ord)
for (i in 1:(ncol(pheno_ord_adni) - 2)) {
  results_ord[[i]][, Phenotype := colnames(pheno_ord_adni)[i + 2]]
}
results_ord <- rbindlist(results_ord)
results_ord[, OR_adni := round(exp(beta), 2)]
results_ord[, P_adni := signif(pval.spa, 2)]
setnames(results_ord, "SNPID", "Variant")

results_bin <- get_results_list(results_files_bin)
for (i in 1:(ncol(pheno_bin) - 2)) {
  phenotype <- colnames(pheno_bin)[i + 2]
  results_bin[[i]][, Phenotype := phenotype]
  results_bin[[i]] <- results_bin[[i]][
    MarkerID %in% results_signif[Phenotype == phenotype, Variant]]
}
results_bin <- rbindlist(results_bin)
setnames(results_bin, c("CHR", "MarkerID"), c("Chromosome", "Variant"))
results_bin[, BETA := -BETA]
results_bin[, OR_adni := round(exp(BETA), 2)]
results_bin[, P_adni := signif(p.value, 2)]

rbind_cols <- c("Phenotype", "Variant", "OR_adni", "P_adni")
results_adni <- rbind(results_bin[, ..rbind_cols],
                     results_ord[, ..rbind_cols])

#-------------------------------------------------
# Merge and format for manuscript
#-------------------------------------------------
results_signif <- merge(results_signif, 
                        results_rosmap, 
                        c("Phenotype", "Variant"),
                        all.x = TRUE)

results_signif <- merge(results_signif, 
                        results_act, 
                        c("Phenotype", "Variant"),
                        all.x = TRUE)

results_signif <- merge(results_signif, 
                        results_adni, 
                        c("Phenotype", "Variant"),
                        all.x = TRUE)

results_signif[, Phenotype := 
              factor(
                Phenotype, 
                labels = c("Arteriolosclerosis", "Atherosclerosis",
                           "Braak NFT Stage", "CAA", "CERAD Score", 
                           "Diffuse Amyloid Plaques", "Gross Infarction",
                           "HS", "LATE-NC", "Lewy Body", "Microinfarct"))]

setnames(results_signif, 
         c("MAF", "OR", "95% CI", "P", "OR_rosmap", "P_rosmap", "OR_act", 
           "P_act", "OR_adni", "P_adni"),
         c("NACC MAF", "NACC OR", "NACC 95% CI", "NACC P-value", 
           "ROSMAP OR", "ROSMAP P-value", "ACT OR", "ACT P-value", 
           "ADNI OR", "ADNI P-value"))
fwrite(results_signif, file = "doc/Supplementary_Table_NACC_GWAS_results.csv",
       quote = F, na = "-")
