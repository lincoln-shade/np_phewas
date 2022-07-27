
library(data.table)
library(rtracklayer)
library(stringi)

genes <- rtracklayer::import('/home/lmsh224/gencode.v40.annotation.gtf.gz')
genes <- as.data.table(genes)
genes[, gene_id_trunc := stri_replace_all_regex(gene_id, '[.].*', '')]

coloc_dir <- 'output/coloc/mega/'
coloc_files <- list.files(coloc_dir)
gwas_phenos <- stri_replace_all_regex(coloc_files, '_.*', '')
coloc_results <- vector(mode = 'list', length = length(coloc_files))

for (i in seq_along(coloc_files)) {
  coloc_results[[i]] <- fread(paste0(coloc_dir, coloc_files[i]), fill = TRUE)
  coloc_results[[i]][, gwas_pheno := ..gwas_phenos[i]]
}

coloc <- rbindlist(coloc_results)

coloc[, gene_id := stri_replace_all_regex(Phenotype, '.*:', '')]
coloc[, gene_id_trunc := stri_replace_all_regex(gene_id, '[.].*', '')]


coloc_genes <- genes[gene_id_trunc %in% coloc$gene_id_trunc
      ][, .(gene_id_trunc, gene_name, gene_type)]
coloc_genes <- coloc_genes[!duplicated(coloc_genes)]
coloc <- merge(coloc, coloc_genes,
                  by = 'gene_id_trunc')

coloc_top <- coloc[PPH4 > 0.9, .(gwas_pheno,gene_name, Tissue, SNP, PPH4, 
                                 beta_GWAS, P_GWAS, beta_QTL, P_QTL)]
setorder(coloc_top, gwas_pheno)
formatC(hs_top$beta_GWAS, digits = 2)

fwrite(coloc_top, file = 'output/coloc/mega_top_coloc_results.csv')
