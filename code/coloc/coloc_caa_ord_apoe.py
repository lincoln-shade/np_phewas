# run in coloc directory
import subprocess
import os
import pandas as pd

pheno_file = "../np_phewas/data/mega/mega_np_ord_dichotomized.pheno"
pheno_data = pd.read_csv(pheno_file)
phenos = ["caa_ord"]

for pheno in phenos:
    cmd = [
        "python", "./coloc_qtl.py",
        "--out_pref", "../np_phewas/output/coloc/mega/" + pheno + "_apoe/",
        "--gwas_sumstats",
        "../np_phewas/output/gwas/mega/polmm/caa_ord_apoe_polmm_results.txt",
        "--max_pval", "1e-5",
        "--pval_col", "pval.spa",
        "--phenotype_type", "cc",
        "--phenotype_col", pheno,
        "--phenotype_file", pheno_file,
        "--phenotype_na_string", "NA",
        "--gtex",
        "--rosmap"
    ]

    subprocess.run(cmd)

