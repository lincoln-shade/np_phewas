#!/usr/bin/python3

import argparse
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True, help = "vcf file path")
parser.add_argument("--geno", help = "maximum propotion of people can be missing alleles")
parser.add_argument("--maf", help = "minimum minor allele frequency") 
parser.add_argument("--out", default = "data/tmp/tmp", help = "output fileset prefix")

args = parser.parse_args()

cmd = ["plink", "--vcf-filter", "--vcf", args.vcf, "--geno", args.geno, "--make-bed", "--biallelic-only", "strict", \
"--set-missing-var-ids", "@:#[b37]:\$1:\$2", "--out", args.out] # "--maf", args.maf,

subprocess.run(cmd)
