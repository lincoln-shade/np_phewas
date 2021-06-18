#!/usr/bin/env python
## updated 1/25/2021

### convert 1000g vcf to plink files
def convert_1kg(i):
    print("- Converting 1000g vcf to PLINK files" + " in chr" + str(i))
    vcf_list = glob.glob("/data_global/1000g/" + args.build + "/" + args.superpop + "/*.vcf.gz")
    chrnum = "chr" + str(i) + ".phase3"
    vcf_file = [s for s in vcf_list if chrnum in s][0]
  
    cmd_convert_1kg = ["plink", "--vcf", str(vcf_file), "--extract", cwd + "/pca_tmp/prune.1000genomes.no_at_cg.in", "--make-bed", "--out", cwd + "/pca_tmp/" + args.superpop + ".chr" + str(i), "--vcf-half-call", "m"]
    subprocess.call(cmd_convert_1kg, stdout=devnull, stderr=devnull)
    
    ### create merge.list for 1000g
    if i > 1:
        f = open(cwd + "/pca_tmp/1000g.merge.list", "a")
        f.write(cwd + "/pca_tmp/" + args.superpop + ".chr" + str(i) + "\n")

## multip threads
def mp(k):
    processes = []
    for j in range(int(args.threads)):
        if int(args.threads)*k + j + 1 == 23:
            break
        p = multiprocessing.Process(target = convert_1kg, args = (int(args.threads)*k + j + 1,))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

### pruning SNPS in the data         
def prunein():
    print("- Pruning SNPs in your data")
    filename = os.path.splitext(ntpath.basename(args.bfile))[0]
    cmd_prunein = ["plink", "--bfile", args.bfile, "--indep-pairwise", "2000", "50", "0.1", "--out", cwd + "/pca_tmp/data"]
       
    subprocess.call(cmd_prunein, stdout=devnull, stderr=devnull)

def no_at_cg():
    print("- Extracting overlapped non-AT-CG SNPs between your data and 1000g")
    
    if args.prune_in:
        prunein = args.prune_in
    else:
        prunein = cwd + "/pca_tmp/data.prune.in"
        
    #cmd_extract = ["cat", "<(cat", abs_path + "/ALL.1000genemes.no-at-cg.uniq.rs)", "<(cat", prunein, "|", "cut", "-f", "4)", "|", "sort", "|", "uniq -d", ">", cwd + "/pca_tmp/prune.1000genomes.no_at_cg.in"]
        
    cmd_extract = "cat <(cat " + abs_path + "/ALL.1000genemes.no-at-cg.uniq.rs) <(cat " + prunein + " | cut -f 4) | sort | uniq -d >" + cwd + "/pca_tmp/prune.1000genomes.no_at_cg.in"
    
    ### List format cannot work
    subprocess.call(cmd_extract, shell = True, executable='/bin/bash')            
            

def plink_no_at_cg(outname):
    print("- Generating pruned PLINK files for your data")
    filename = os.path.splitext(ntpath.basename(args.bfile))[0]
    cmd_prunedfile = ["plink", "--bfile", args.bfile, "--extract", cwd + "/pca_tmp/prune.1000genomes.no_at_cg.in", "--make-bed", "--out", cwd + "/pca_tmp/pruned." + outname]
    subprocess.call(cmd_prunedfile, stdout=devnull, stderr=devnull)
        
def merge_plink_1000g(outname):
    print("- Merging 1000g and your data")
    cmd_merge_plink_1000g = ["plink", "--bfile", cwd + "/pca_tmp/pruned." + outname, "--bmerge", cwd + "/pca_tmp/" + args.superpop, "--make-bed", "--out", cwd + "/pca_tmp/pruned." + outname + "." + args.superpop]
    subprocess.call(cmd_merge_plink_1000g, stdout=devnull, stderr=devnull)

    if os.path.exists(cwd + "/pca_tmp/pruned." + outname + "." + args.superpop + "-merge.missnp"):
        print("Correcting strand error")
        cmd_re1 = ["plink", "--bfile", cwd + "/pca_tmp/pruned." + outname, "--exclude", cwd + "/pca_tmp/pruned." + outname + "." + args.superpop + "-merge.missnp", "--make-bed", "--out", cwd + "/pca_tmp/pruned." + outname]
        cmd_re2 = ["plink", "--bfile", cwd + "/pca_tmp/" + args.superpop, "--exclude", cwd + "/pca_tmp/pruned." + outname + "." + args.superpop + "-merge.missnp", "--make-bed", "--out", cwd + "/pca_tmp/" + args.superpop]
        subprocess.call(cmd_re1, stdout=devnull, stderr=devnull)
        subprocess.call(cmd_re2, stdout=devnull, stderr=devnull)
        subprocess.call(cmd_merge_plink_1000g, stdout=devnull, stderr=devnull)
    

def output_pcs(outname):    
    evec = pd.read_csv(cwd + "/pca_tmp/pruned." + outname + "." + args.superpop + ".pca.eigenvec", delim_whitespace=True, header = None)
    fam = pd.read_csv(cwd + "/pca_tmp/pruned." + outname + "." + args.superpop + ".fam", delim_whitespace=True, header = None, na_values = "NA")
    pop = pd.read_csv("/data_global/1000g/integrated_call_samples_v3.20130502.ALL.panel", delim_whitespace=True)
    
    pop = pop.iloc[:,[0,2]]
    fam = fam.iloc[:,[0,1]]
    pop = pop.rename(columns={"sample": "IID", "super_pop": "Superpop"})
    fam = fam.rename(columns={0: "FID", 1: "IID"})
    evec.columns = ['FID', 'IID'] + ["V" + str(s) for s in list(range(1,21,1))]
    
    dat = pd.merge(fam, pop, on='IID', how='left')
    dat = pd.merge(dat, evec, on=["FID", 'IID'], how='left')
    dat.to_csv(cwd + "/" + outname + ".pcs.csv", index = False)


def umap_pca(outname):
    dat = pd.read_csv(cwd + "/" + outname + ".pcs.csv", delim_whitespace=True)
    dat_pc = dat.iloc[:,2:22]
    tmp = dat_pc.values
    reducer = umap.UMAP(random_state = 7, n_components=2)  
    embedding = reducer.fit_transform(tmp)   
    embedding.shape
    dat_pc = pd.concat([dat.iloc[:,[0,1]], pd.DataFrame(embedding[:,:2]), dat.iloc[:,[22,23]]], axis=1)
    dat_pc = dat_pc.rename(columns={0: "V1", 1: "V2"})
    dat_pc.to_csv(cwd + "/" + outname + ".embedding.csv", index = False)
    
    
if __name__ == '__main__':    
    
    import argparse
    import subprocess
    import os
    import multiprocessing
    import ntpath
    import math   
    import glob
    import pandas as pd
    import numpy as np
    #import umap
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--bfile", help = "plink binary file")
    parser.add_argument("--superpop", default = "ALL", help = "specify super population")
    parser.add_argument("--threads", default = '1', help = "the number of thread to use")
    parser.add_argument("--build", default = 'hg19', help = "hg18 or hg19 or hg38")  
    parser.add_argument("--out",  default = 'PLINK', help = "output file name")  
    parser.add_argument("--prune-in",  help = "prune in file if you have") 
    parser.add_argument("--umap", action = "store_true", help = "wheter umap is performed") 
    args = parser.parse_args()
    
    ## path to this file
    abs_path = os.path.dirname(os.path.abspath(__file__))

   ## output path and set the file name
    try:        
        cwd = os.path.dirname(os.path.abspath(args.out))
        plinkfile = ntpath.basename(args.out)
    except:
        cwd = os.getcwd()
        plinkfile = ntpath.basename(args.bfile)
 
    # make pca_tmp directory for working
    if not os.path.exists(cwd + "/pca_tmp/"):
        os.mkdir(cwd + "/pca_tmp/")    
        
    ## 
    if os.path.exists(cwd + "/pca_tmp/1000g.merge.list"):
        os.remove(cwd + "/pca_tmp/1000g.merge.list")
    
    if os.path.exists(cwd + "/pca_tmp/data.merge.list"):
        os.remove(cwd + "/pca_tmp/data.merge.list")         
    
    devnull = open(os.devnull, 'wb')

    ### run prunein
    prunein()
    
    ### run no_at_cg
    no_at_cg()
    
    ## run plink_no_at_cg()
    plink_no_at_cg(plinkfile)

    ### run mp <--- convert_1kg with multi threads
    rep = int(math.ceil(float(22)/int(args.threads)))
    for t in range(rep):
        mp(t)

    ##  merge All 1000 genomes files
    cmd_merge_1000g = ["plink", "--bfile", cwd + "/pca_tmp/" + args.superpop + ".chr1", "--merge-list", cwd + "/pca_tmp/1000g.merge.list", "--make-bed", "--out", cwd + "/pca_tmp/" + args.superpop]
    subprocess.call(cmd_merge_1000g, stdout=devnull, stderr=devnull)   

    ## run merge_plink_1000g()    
    merge_plink_1000g(plinkfile)    
       
    cmd_pedind = ["Rscript", "--slave", "--vanilla", abs_path + "/pedind.r", args.superpop, cwd + "/pca_tmp/pruned." + plinkfile + "." + args.superpop]
    subprocess.call(cmd_pedind, stdout=devnull, stderr=devnull)

    cmd_pca = ["plink", "--bfile", cwd + "/pca_tmp/pruned." + plinkfile + "." + args.superpop, "--pca", "--out", cwd + "/pca_tmp/pruned." + plinkfile + "." + args.superpop + ".pca", "--allow-no-sex", "--threads", str(args.threads)]
    print("- Calculating PCs")
    subprocess.call(cmd_pca, stdout=devnull, stderr=devnull)

    output_pcs(plinkfile)

    ## create a scatter plot of 1st and 2nd PCs
    print("- Creating PC plot")
    cmd_plot = ["Rscript", "--slave", "--vanilla", abs_path + "/plot-pca-results.r", cwd + "/" + plinkfile + ".pcs.csv"]
    subprocess.call(cmd_plot)
    
    if args.umap:
        umap_pca(plinkfile)
        cmd_plot = ["Rscript", "--slave", "--vanilla", abs_path + "/plot-pca-results.r", cwd + "/" + plinkfile + ".embedding.csv"]
        subprocess.call(cmd_plot, stdout=devnull, stderr=devnull)
    
    
    print("(~^.^)~ Done ~(^.^~)")
    

