library(data.table)
library(rtracklayer)
library(stringi)
# colocalization results
pheno <- fread("shared/nacc_rosmap_act_np.pheno")
phenos <- colnames(pheno)[3:ncol(pheno)]
coloc_list <- vector(mode = "list", length = length(phenos))
names(coloc_list) <- phenos
for (i in 1:length(phenos)) {
  p <- phenos[i]
  coloc_path <- paste0("output/coloc/meta/", p, "/coloc_results.txt")
  coloc_list[[i]] <- fread(coloc_path)
}

coloc_list[[length(coloc_list) + 1]] = fread("output/coloc/meta/caa_apoe/coloc_results.txt")
names(coloc_list)[[length(coloc_list)]] = "caa_apoe_adj"
coloc <- rbindlist(coloc_list)
coloc[, gene_id := stri_replace_all_regex(QTL_Phenotype, '.*:', '')]
coloc[, gene_id := stri_replace_all_regex(gene_id, '\\..*', '')]

# gene names
genes <- rtracklayer::import('/home/lmsh224/downloads/gencode.v40.annotation.gtf.gz')
genes <- as.data.table(genes)
genes[, gene_id_trunc := stri_replace_all_regex(gene_id, '\\..*', '')]

coloc_genes <- genes[gene_id_trunc %in% coloc$gene_id
][, .(gene_id_trunc, gene_name, gene_type)]
coloc_genes <- coloc_genes[!duplicated(coloc_genes)]
setnames(coloc_genes, "gene_id_trunc", "gene_id")

coloc <- merge(coloc, coloc_genes, by = 'gene_id', all.x = TRUE)
coloc[is.na(gene_name), gene_name := QTL_Phenotype]
setnames(coloc, "gene_name", "QTL_Target")

fwrite(coloc, "output/coloc/meta/coloc_results.txt")
fwrite(coloc[PPH4 > 0.8], "output/coloc/meta/coloc_results_pph4_80.csv")
