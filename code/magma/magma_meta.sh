
for pheno in diffuse_abeta
    do
    
    n=$(./code/magma/count_non_missing_values.R shared/nacc_rosmap_act_np.pheno $pheno)
    echo $n
    
    magma_input=tmp/"$pheno"_magma_input.csv
    ./code/magma/extract_markername_pvalue.R \
        -i output/gwas/metal/results/"$pheno"1.csv \
        -o "$magma_input"
    
    head "$magma_input" 
    
    magma \
        --bfile data/adc/adc_np \
        --gene-annot data/adc/snp_gene_annot_hg38.genes.annot \
        --pval "$magma_input" \
        N="$n" \
        --out output/magma/"$pheno"
    
    done

