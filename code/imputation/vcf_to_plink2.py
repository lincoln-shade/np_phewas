#!/usr/bin/python3

def vcf_to_plink2(vcf):
    sys.stdout = sys.__stdout__
    
    cmd = ['plink2', '--vcf', vcf, '--geno', args.geno, \
    '--maf', args.maf, '--make-pgen', "--out", args.out, '--vcf-half-call', 'm']
    
    subprocess.run(cmd)

if __name__ == '__main__':

    import argparse
    import subprocess
    import sys
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', 
                        required = True, 
                        help = "vcf file path")
    parser.add_argument('--geno', 
                        default = '0.2', 
                        help = "max fraction of people can be missing alleles")
    parser.add_argument('--maf', 
                        default = '0.001', 
                        help = "minimum minor allele frequency") 
    parser.add_argument('--out', 
                        default = 'data/tmp/tmp', 
                        help = "output fileset prefix")

    args = parser.parse_args()
    
    cwd = os.getcwd()

    vcf_file = args.vcf
    
    vcf_to_plink2(vcf_file)    
