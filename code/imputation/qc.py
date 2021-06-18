#!/usr/bin/env python
## updated 2021-06-13 by lincoln-shade


# def ibs(filename):
# 
#     ## read .genome and randomly obtain one of pair
#     ibs_df = pd.read_csv(cwd + rel_tmp_path + filename + ".genome", sep ="\s*", engine='python')
#     ibs_df = ibs_df[["FID1", "IID1", "FID2", "IID2"]]
#     remove_ids = []
#     for index, dat in ibs_df.iterrows():
#         rand_num = random.choice([0,2]) ## random value 0 or 2
#         remove_ids.append([dat[rand_num], dat[(rand_num + 1)]])
#     remove_ids = pd.DataFrame(remove_ids).drop_duplicates()
#     remove_ids.to_csv(cwd + rel_tmp_path + filename + ".genome.removed.ids", header = False, index = False, sep=" ")
#     if os.path.exists(cwd + rel_aux_path + filename + "_genome_removed.ids"): 
#         remove_ids.to_csv(cwd + rel_aux_path + filename + "_genome_removed.ids", mode = "a", header = False, index = False, sep=" ")
#     else:
#         remove_ids.to_csv(cwd + rel_aux_path + filename + "_genome_removed.ids", mode = "w", header = False, index = False, sep=" ")
    
        

def qc(filename):
    sys.stdout = sys.__stdout__
    ## 1. remove samples that are missing more than 5%
    print("Setp 1: Removing subjects with missing rate > " + args.mind)
    cmd_step1 = ["plink", "--bfile", args.bfile, "--mind", args.mind, "--make-bed", "--out", cwd + rel_tmp_path + filename + "_mind"]
    if args.keep:
        cmd_step1 = cmd_step1 + ["--keep", args.keep]
    
    subprocess.call(cmd_step1, stdout=devnull, stderr=devnull)

        
    ## output removed sample IDs
    cmd_step1_auxiliary = ["cat", args.bfile + ".fam", cwd + rel_tmp_path + filename + "_mind.fam|sort|uniq", "-d>", cwd + rel_aux_path + filename + "_mind_removed.ids"]

    subprocess.call(cmd_step1_auxiliary, stdout=devnull, stderr=devnull)


    ## 2. remove snp with MAF < args.maf (if present)
    if args.maf:
        print("Step 2: Removing SNPs with MAF < " + args.maf)
        cmd_step2 = ["plink", "--bfile", cwd + rel_tmp_path + filename + "_mind", "--maf", args.maf, "--make-bed", "--out", cwd + rel_tmp_path + filename + "_maf"]
        subprocess.call(cmd_step2, stdout=devnull, stderr=devnull)
    
        ## output removed SNP IDs
        cmd_step2_auxiliary = ["cat", args.bfile + ".bim", cwd + rel_tmp_path + filename + "_maf.bim|sort|uniq", "-d>", cwd + rel_aux_path + filename + "_maf_removed.snps"]

        subprocess.call(cmd_step2_auxiliary, stdout=devnull, stderr=devnull)
 
    ## 3. remove snps that are missing more than 5%
    print("Step 3: Removing SNPs with missing rate > " + args.geno)
    if args.maf:
        suffix = "_maf"
    else:
        suffix = "_mind"
    cmd_step3 = ["plink", "--bfile", cwd + rel_tmp_path + filename + suffix, "--geno", args.geno, "--make-bed", "--out", cwd + rel_tmp_path + filename + "_geno"]
    subprocess.call(cmd_step3, stdout=devnull, stderr=devnull)
    
    ## output removed SNP IDs
    cmd_step3_auxiliary = ["cat", cwd + rel_tmp_path + filename + suffix + ".bim", cwd + rel_tmp_path + filename + "_geno.bim|sort|uniq", "-d>", cwd + rel_aux_path + filename + "_geno_removed.snps"]
    subprocess.call(cmd_step3_auxiliary, stdout=devnull, stderr=devnull)

    # ## 4. pruning SNPs based on LD
    # print("Step 4: Pruning SNPs")
    # cmd_step4 = ["plink", "--bfile", cwd + rel_tmp_path + filename + "_geno", "--indep-pairwise", "1500", "50", "0.2", "--out", cwd + rel_tmp_path + filename]
    # subprocess.call(cmd_step4, stdout=devnull, stderr=devnull)
   
    
    # ## 5. remove one of pairs with PHI_HAT > 0.185
    # print("Step 5: Removing one of pair with PHI_HAT > " + args.ibd)
    # 
    # ## delete _genome_removed.ids if it exists
    # if os.path.exists(cwd + rel_aux_path + filename + "_genome_removed.ids"):
    #     os.remove(cwd + rel_aux_path + filename + "_genome_removed.ids")
    # 
    # ## copy _geno to _ibs
    # copyfile(cwd + rel_tmp_path + filename + "_geno.bim", cwd + rel_tmp_path + filename + "_ibs.bim")
    # copyfile(cwd + rel_tmp_path + filename + "_geno.bed", cwd + rel_tmp_path + filename + "_ibs.bed")
    # copyfile(cwd + rel_tmp_path + filename + "_geno.fam", cwd + rel_tmp_path + filename + "_ibs.fam")
    # 
    # ## calculate ibs (PHI_HAT) 
    # cmd_step5_1 = ["plink", "--bfile", cwd + rel_tmp_path + filename + "_ibs", "--extract", cwd + rel_tmp_path + filename + ".prune.in", "--genome", "--min", args.ibd, "--out",  cwd + rel_tmp_path + filename]
    #     
    # subprocess.call(cmd_step5_1, stdout=devnull, stderr=devnull)    
    
    ## create cwd + rel_tmp_path + filename + ".genome.removed.ids"
    # ibs(filename)
    
    # while os.stat(cwd + rel_tmp_path + filename + ".genome.removed.ids").st_size > 0: ## loop until all are independent
    #     
    #     ## run plink to remove subjects with PHI_HAT > 0.185
    #     cmd_step5_2 = ["plink", "--bfile", cwd + rel_tmp_path + filename + "_ibs", "--remove", cwd + rel_tmp_path + filename + ".genome.removed.ids", "--make-bed", "--out",  cwd + rel_tmp_path + filename + "_ibs"]
    # 
    #     subprocess.call(cmd_step5_2, stdout=devnull, stderr=devnull)
    #     subprocess.call(cmd_step5_1, stdout=devnull, stderr=devnull)
    # 
    #     ibs(filename)    
    
    ## 6. remove samples based inbreeding coefficients
    print("Step 6: Removing subjects based on inbreeding coefficients")
    ## calculate inbreeding coefficients
    cmd_step6_1 = ["plink", "--bfile", cwd + rel_tmp_path + filename + "_geno", "--het", "--out", cwd + rel_tmp_path + filename]
    subprocess.call(cmd_step6_1, stdout=devnull, stderr=devnull)
     
    ## import .het and extract sample with F > mean F + 3SD F or F < mean F - 3SD F in het.r > remove.het
    het_df = pd.read_csv(cwd + rel_tmp_path + filename + ".het", engine='python', delim_whitespace = True)
    het_lower = het_df["F"].mean() - 3*het_df["F"].std()
    het_upper = het_df["F"].mean() + 3*het_df["F"].std()

    het_df = het_df.loc[(het_df["F"] < het_lower) & (het_df["F"] > het_upper)]
    if het_df.empty==False:
        het_df.to_csv(cwd + rel_aux_path + filename + "_het_removed.ids", sep = " ", index = False, header = False)
    
        ## remove samples
        cmd_step6_2 =  ["plink", "--bfile", cwd + rel_tmp_path + filename + "_geno", "--remove",  cwd + rel_aux_path + filename + "_het_removed.ids", "--make-bed", "--out", cwd + rel_tmp_path + filename + ".qced"]
        subprocess.call(cmd_step6_2, stdout=devnull, stderr=devnull)
    else:
        copyfile(cwd + rel_tmp_path + filename + "_geno.bim", cwd + rel_tmp_path + filename + ".qced.bim")
        copyfile(cwd + rel_tmp_path + filename + "_geno.bed", cwd + rel_tmp_path + filename + ".qced.bed")
        copyfile(cwd + rel_tmp_path + filename + "_geno.fam", cwd + rel_tmp_path + filename + ".qced.fam")
        
    print("Done! o(^_^)/~~")
    
if __name__ == '__main__':

    import argparse
    import subprocess
    import os
    import multiprocessing
    import math
    import pandas as pd
    import random
    import sys
    from shutil import copyfile
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--bfile", required=True, help = "plink binary file name") ## plink file path
    parser.add_argument("--mind", nargs = "?", default = "0.2", help = "missing rate to remove individuals. Default is 0.05") ## e.g., 0.05
    parser.add_argument("--maf", help = "MAF to remove markers") ## e.g., 0.05
    parser.add_argument("--geno", nargs = "?", default = "0.2", help = "missing rate to remove markers. Default is 0.05") ## e.g., 0.05
    # parser.add_argument("--ibd", help = "exclusion criteria to remove one of pairs. Default is 0.185") ## e.g., 0.25
    parser.add_argument("--keep", help = "ID file with FID and IID to extract a subset of individuals")
    parser.add_argument("--pheno", help = "phenotype file with FID, IID, binary phenotype (0 = unaffected, 1 = affected) if case-control study")
    parser.add_argument("--pheno-name", help = "phenotype name in case-control study if --pheno is specified and phenotype is not on the third column")
    parser.add_argument("--out", help = "output QCed plink file name")

    args = parser.parse_args()

    ## path to this file
    abs_path = os.path.dirname(os.path.abspath(__file__))

   ## output path and set the file name
    try:        
        # cwd = os.path.dirname(os.path.abspath(args.out))
        cwd = os.getcwd()
        rel_tmp_path = "/data/tmp/"
        rel_aux_path = "/data/auxiliary"
        plinkfile = os.path.basename(args.out)
    except:
        plinkfile = os.path.basename(args.bfile)
    
    
    
    # make tmp directory for working
    if not os.path.exists(cwd + rel_tmp_path):
        os.mkdir(cwd + rel_tmp_path)    
    
     # make auxiliary directory
    if not os.path.exists(cwd + rel_aux_path):
        os.mkdir(cwd + rel_aux_path)    
        
    devnull = open(os.devnull, 'wb')    
    qc(plinkfile)    
