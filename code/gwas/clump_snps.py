#!/usr/bin/python3

import sys
import argparse
import os
import subprocess

def clump_snps(results, 
               outfile,
               bfile,
               p1 = None, 
               r2 = None, 
               snp_field = None,
               p_field = None):
    
    cmd = ["plink", "--clump", results, "--out", outfile, "--bfile", bfile]
    
    if p1:
        cmd = cmd + ["--clump-p1", p1]
    
    if r2:
        cmd = cmd + ["--clump-r2", r2]
    
    if snp_field:
        cmd = cmd + ["--clump-snp-field", snp_field]
    
    if p_field:
        cmd = cmd + ["--clump-field", p_field]
    
    subprocess.run(cmd)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--results', 
                        help = "file with SNPs and p-values")
    parser.add_argument('-o', '--out', 
                        help = "prefix for clump output file")
    parser.add_argument('-b', '--bfile', 
                        help = "prefix of PLINK binary fileset")
    parser.add_argument('-p', '--p1', default = '0.00001',
                        help = "maximum p-value for lead clump SNP")
    parser.add_argument('-r2', '--r2', default = "0.05",
                        help = "minimum LD r2 for clumping")
    parser.add_argument('--snp_field', default = 'SNP', 
                        help = "SNP ID column name")
    parser.add_argument('--p_field', default = 'P', 
                        help = "p-value column name")
    args = parser.parse_args()
    
    clump_snps(args.results, 
               args.out,
               args.bfile,
               args.p1, 
               args.r2, 
               args.snp_field, 
               args.p_field)
    
    
    
    
    
