phenos=$(cat data/mega/mega_np.pheno | \
         awk 'NR==1{for (i=3; i<NF; i++) printf $i "\n"; print $NF}')
         
for pheno in $phenos
do
  Rscript --vanilla code/magma/02_subset_pheno.R \
    -p data/mega/mega_np.pheno \
    -o tmp/"$pheno".tmp \
    -r data/mega/related_rm/"$pheno".remove \
    --phenotype "$pheno"
done

nprocs=29
for pheno in $phenos
do
  seq 1 $nprocs | xargs -n 1 -P $nprocs -I % \
    magma \
    --batch % $nprocs \
    --bfile data/mega/mega_np \
    --pheno file=tmp/"$pheno".tmp \
    --covar file=data/mega/mega_np.covar \
    --gene-annot data/mega/magma/snp_gene_annot_hg38.genes.annot \
    --out tmp/"$pheno"_magma

  magma --merge tmp/"$pheno"_magma --out output/magma/"$pheno"
done

# ordinal phenotypes
phenos=$(cat data/mega/mega_np_ord.pheno | \
         awk 'NR==1{for (i=3; i<NF; i++) printf $i "\n"; print $NF}')
         
for pheno in $phenos
do
  Rscript --vanilla code/magma/02_subset_pheno.R \
    -p data/mega/mega_np_ord.pheno \
    -o tmp/"$pheno".tmp \
    -r data/mega/related_rm/"$pheno".remove \
    --phenotype "$pheno"
done

for pheno in $phenos
do
  seq 1 $nprocs | xargs -n 1 -P $nprocs -I % \
    magma \
    --batch % $nprocs \
    --bfile data/mega/mega_np \
    --pheno file=tmp/"$pheno".tmp \
    --covar file=data/mega/mega_np.covar \
    --gene-annot data/mega/magma/snp_gene_annot_hg38.genes.annot \
    --genes-only \
    --out tmp/"$pheno"_magma

  magma --merge tmp/"$pheno"_magma --out output/magma/"$pheno"
done