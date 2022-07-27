
# annotate variants to genes for MAGMA gene-based analysis
# build: hg38
# gene window: 100 Mb

cat raw_data/NCBI38.gene.loc | \
  awk '{print $6, $2, $3, $4, $5, $1}' \
  > tmp/hg38.gene.loc

magma \
  --snp-loc data/mega/mega_np.bim \
  --gene-loc tmp/hg38.gene.loc \
  --annotate window=100 \
  --out data/mega/snp_gene_annot_hg38