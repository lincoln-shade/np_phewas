#!/bin/python3

import sys
import argparse
import os
import subprocess

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bfile')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-m', '--misspheno', nargs='?', default = '-1')
    parser.add_argument('-r', '--remove')
    parser.add_argument('-p', '--pheno')
    parser.add_argument('-c', '--covar')
    parser.add_argument('--chrom')
    parser.add_argument('--pheno_name')
    args = parser.parse_args()
    
    cmd = ["plink", "--bfile", args.bfile, "--out", args.outfile, 
    "--missing-phenotype", args.misspheno, "--pheno", args.pheno, 
    "--covar", args.covar, "--1", "--logistic", "hide-covar",
    "--ci", "0.95", "--allow-no-sex"]
    
    if args.chrom:
        cmd = cmd + ["--chr", args.chrom]
    
    if args.pheno_name:
        cmd = cmd + ["--pheno-name", args.pheno_name]
    
    if args.remove:
        cmd = cmd + ["--remove", args.remove]
    
    subprocess.run(cmd)
