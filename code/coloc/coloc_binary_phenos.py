
# run in coloc directory
import subprocess

for pheno in ['grossinf', 'microinf']:
    cmd = [
        "python", "./coloc_qtl.py",
        "--out_pref", "../np_phewas/output/coloc/mega/" + pheno + "/",
        "--gwas_sumstats", 
        "../np_phewas/output/gwas/mega/saige/" + pheno + "_saige_results.txt",
        "--max_pval", "1e-5", 
        "--pval_col", "p.value",
        "--chrom_col", "CHR",
        "--snp_col", "MarkerID",
        "--maf_col", "AF_Allele2",
        "--phenotype_type", "cc",
        "--phenotype_col", pheno,
        "--phenotype_file", "../np_phewas/data/mega/mega_np.pheno",
        "--phenotype_na_string", "-1", 
        "--gtex", 
        "--rosmap"
    ]
    
    subprocess.run(cmd)
