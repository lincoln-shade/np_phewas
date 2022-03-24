#!/usr/bin/python3

import sys
import argparse
import os
import subprocess

def run_plink_regression(pfile, outfile, misspheno, pheno, covar, 
                         chrom = None, remove = None, pheno_name = None,
                         extract = None, keep = None):
    
    cmd = ["plink2", "--pfile", pfile, "--out", outfile, 
    "--input-missing-phenotype", args.misspheno, "--pheno", pheno, 
    "--covar", covar, "--1", "--mac", "20", "--glm", "hide-covar",
    "no-firth", "--ci", "0.95", "--covar-variance-standardize"]
    
    if chrom:
        cmd = cmd + ["--chr", chrom]
    
    if pheno_name:
        cmd = cmd + ["--pheno-name", pheno_name]
    
    if remove:
        cmd = cmd + ["--remove", remove]
    
    if extract:
        cmd = cmd + ["--extract", extract]
    
    if keep:
        cmd = cmd + ["--keep", keep]
    
    subprocess.run(cmd)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--pfile')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-m', '--misspheno', nargs='?', default = '-1')
    parser.add_argument('-r', '--remove')
    parser.add_argument('-p', '--pheno')
    parser.add_argument('-c', '--covar')
    parser.add_argument('--chrom')
    parser.add_argument('--pheno_name')
    parser.add_argument('--extract')
    parser.add_argument('--keep')
    args = parser.parse_args()
    
    if args.chrom:
        chrom = args.chrom
    else:
        chrom = None
    
    if args.remove:
        remove = args.remove
    else:
        remove = None
     
    if args.pheno_name:
        pheno_name = args.pheno_name
    else:
        pheno_name = None
    
    if args.extract:
        extract = args.extract
    else:
        extract = None
    
    if args.keep:
        keep = args.keep
    else:
        keep = None
    
    run_plink_regression(args.pfile, args.outfile, args.misspheno, args.pheno,
                         args.covar, chrom, remove, pheno_name, extract, keep)
