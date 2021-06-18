#!/usr/bin/env python
## updated 6/10/2020

import argparse
import subprocess
import os
import multiprocessing
import math

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--bfile", help = "plink binary file name") ## plink file path
    parser.add_argument("--keep", help = "ID file with FID and IID to keep a subset of individuals")
    parser.add_argument("--remove", help = "ID file with FID and IID to remove a subset of individuals")


    args = parser.parse_args()

    ## current path
    cwd = os.getcwd()
    
    # ## path to this file
    # abs_path = os.path.dirname(os.path.abspath(__file__))

    ## make tmp directory for working under the current directory
    rel_path = "data/tmp/"
    if not os.path.exists(cwd + "/" + rel_path):
        os.mkdir(cwd + "/" + rel_path)
    if not os.path.exists(cwd + "/data/vcf/"):
        os.mkdir(cwd + "/data/vcf/")

    ## extract filename only from full path of plink file
    plinkfile = rel_path + os.path.basename(args.bfile)
    
    ## change working directory to tmp/
    # os.chdir(cwd + "/tmp/")

    ## https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/
    ## 1. keep or/and remove samples
    cmd_step1 = "plink --bfile " + args.bfile + " --make-bed --out " + plinkfile + "_extract"
    if args.keep:
        keep = " --keep " + args.keep
        cmd_step1 = cmd_step1 + keep
    if args.remove:
        remove = " --remove " + args.remove
        cmd_step1 = cmd_step1 + remove
    subprocess.call(cmd_step1, shell = True)

    ## 2. Create a frequency file
    cmd_step2 = "plink --freq --bfile " + plinkfile + "_extract --out " + plinkfile + "_extract"
    subprocess.call(cmd_step2, shell = True)

    ## 3. Execute script
    cmd_perl = "perl " + "code/imputation/HRC-1000G-check-bim.pl -b " + plinkfile + "_extract.bim -f " + plinkfile + "_extract.frq -r " + "raw_data/PASS.Variantsbravo-dbsnp-all.tab.gz -h"
    subprocess.call(cmd_perl, shell = True)

    ## 4. run shell script
    ## (perl script writes PLINK commands to this script)
    cmd_sh = "sh " + rel_path + "Run-plink.sh"
    subprocess.call(cmd_sh, shell = True)

    ## 5. convert plink to vcf
    for i in range(22):
        cmd_step5 = "plink --bfile " + plinkfile + "_extract-updated-chr" + str(i+1) + " --recode vcf --out " + plinkfile + "_extract-updated-chr" + str(i+1)
        subprocess.call(cmd_step5, shell = True)
    
    ## 6. create vcf.gz
    for i in range(22):
        cmd_step6 = "vcf-sort " + plinkfile + "_extract-updated-chr" + str(i+1) + ".vcf | bgzip -c > data/vcf/" + os.path.basename(args.bfile) + "_extract-updated-chr" + str(i+1) + ".vcf.gz"
        subprocess.call(cmd_step6, shell = True)

    
    
    
    
