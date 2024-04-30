# run in coloc directory
import subprocess
import os
import pandas as pd

pheno_file = "../np_phewas/shared/nacc_rosmap_act_np_dichotomized.pheno"
pheno_data = pd.read_csv(pheno_file)
phenos = ["caa"]

for pheno in phenos:
    cmd = [
        "python", "coloc_qtl.py",
        "--out_pref", "../np_phewas/output/coloc/meta/" + pheno + "_apoe/",
        "--gwas_sumstats", 
        "../np_phewas/output/gwas/metal/results/" + pheno + "_apoe1.csv",
        "--max_pval", "1e-5",
        "--pval_col", "P.value",
        "--maf_col", "Freq1",
        "--chrom_col", "chr",
        "--snp_col", "MarkerName",
        "--beta_col", "Effect",
        "--betasd_col", "StdErr",
        "--phenotype_type", "cc",
        "--phenotype_col", pheno,
        "--phenotype_file", pheno_file,
        "--phenotype_na_string", "NA", 
        "--gtex", 
        "--rosmap"
    ]

    subprocess.run(cmd)

