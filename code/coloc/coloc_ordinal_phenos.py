# run in coloc directory
import subprocess
import os
import pandas as pd

pheno_file = "../np_phewas/data/mega/mega_np_ord_dichotomized.pheno"
pheno_data = pd.read_csv(pheno_file)
phenos = pheno_data.columns[2:]

for pheno in phenos:
    cmd = [
        "python", "./src/coloc_qtl.py",
        "--out_pref", "../np_phewas/output/coloc/mega/" + pheno + "/",
        "--gwas_sumstats", 
        "../np_phewas/output/gwas/mega/polmm/" + pheno + "_polmm_results.txt",
        "--max_pval", "1e-5",
        "--pval_col", "pval.spa",
        "--phenotype_type", "cc",
        "--phenotype_col", pheno,
        "--phenotype_file", pheno_file,
        "--phenotype_na_string", "NA", 
        "--gtex", 
        "--rosmap"
    ]
    
    subprocess.Popen(cmd)
