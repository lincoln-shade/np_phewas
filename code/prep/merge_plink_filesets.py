
#============================
# merge two plink filesets
#============================

def tmp_file(bfile, suffix=''):
    return tmp_path + "/" + os.path.basename(bfile) + suffix

def strip_a1a2(bfile, suffix):
    # strips the A1:A2 suffixes from variant rsIDs
    out_path = tmp_file(bfile, suffix)
    with open(out_path, "w") as file_out:
        cat_proc = subprocess.Popen(["cat", bfile + ".bim"], 
                                    stdout=subprocess.PIPE)
        awk_proc = subprocess.Popen(["awk", 
                                     "{print gensub(/:[ACTG]*:[ACTG]*/, \"\", 1)}"],
                                    stdin=cat_proc.stdout,
                                    stdout=file_out)
        awk_proc.wait()

def remove_dup_rsids(bfile, suffix, suffix1):
    bim_file = tmp_file(bfile, suffix)
    out_path = tmp_file(bfile, "_dup_rsids.tmp")
    with open(out_path, "w") as file_out:
        cat_proc = subprocess.Popen(["cat", bim_file], 
                                    stdout=subprocess.PIPE)
        cut_proc = subprocess.Popen(['awk', '{print $2}'], 
                                    stdin=cat_proc.stdout,
                                    stdout=subprocess.PIPE)
        sort_proc = subprocess.Popen(['sort'], 
                                    stdin=cut_proc.stdout,
                                    stdout=subprocess.PIPE)
        uniq_proc = subprocess.Popen(['uniq', '-d'], 
                                    stdin=sort_proc.stdout,
                                    stdout=file_out)
        uniq_proc.wait()
    
    cmd = ['plink', '--bfile', bfile, '--bim', bim_file, '--make-bed',
           '--out', tmp_file(bfile, suffix1)]

    if os.stat(out_path).st_size != 0:
        cmd += ['--exclude', out_path]

    subprocess.run(cmd)

def print_variants(bfile, suffix, suffix1):
    with open(tmp_file(bfile, suffix)) as in_file:
        with open(tmp_file(bfile, suffix1), 'w') as out_file:
            subprocess.run(['awk', '{print$2}'],
                           stdin=in_file,
                           stdout=out_file)

def extract_variants(bfile, variants, suffix):
    cmd = ['plink', '--bfile', bfile, '--extract', variants, 
           '--make-bed', '--out', tmp_file(bfile, suffix)]
    subprocess.run(cmd)
    
def sync_maps(bfile, bfile1, suffix):
    subprocess.run(['plink', '--bfile', bfile1, '--update-map', 
                    bfile + '.bim', '4', '2', '--make-bed', 
                    '--out', bfile1 + suffix])

def sync_alleles(bfile, bfile1, suffix):
    # 1. set bfile1 alleles to bfile alleles
    subprocess.run(['plink', '--bfile', bfile1, '--a1-allele', bfile + '.bim',
                    '5', '2', '--make-bed', '--out', bfile1 + '_upmap'])
    
    # 2. try to resolve strand issues
    with open(bfile + '_vars.tmp', 'w') as out_file:
        subprocess.run(['awk', '{print $2, $5, $6}', bfile + '.bim'],
                       stdout=out_file)
    
    with open(bfile1 + '_vars.tmp', 'w') as out_file:
        subprocess.run(['awk', '{print $2, $5, $6}', bfile1 + '_upmap.bim'],
                       stdout=out_file)
    
    with open(bfile1 + '_flip_list.tmp', 'w') as out_file: 
        sort_proc = subprocess.Popen(['sort', bfile + '_vars.tmp', 
                                      bfile1 + '_vars.tmp'], 
                                     stdout=subprocess.PIPE)
    
        uniq_proc = subprocess.Popen(['uniq', '-u'],
                                     stdin=sort_proc.stdout,
                                     stdout=subprocess.PIPE)
        
        awk_proc = subprocess.Popen(['awk', '{print $1}'],
                                    stdin=uniq_proc.stdout,
                                    stdout=subprocess.PIPE)
        
        sort_proc1 = subprocess.Popen(['sort', '-u'],
                                      stdin=awk_proc.stdout,
                                      stdout=out_file)
        
        sort_proc1.wait()
    
    subprocess.run(['plink', '--bfile', bfile1 + '_upmap', 
                    '--flip', bfile1 + '_flip_list.tmp',
                    '--a1-allele', bfile + '.bim', '5', '2', '--make-bed',
                    '--out', bfile1 + '_flipped'])
    
    with open(bfile1 + '_flipped_vars.tmp', 'w') as out_file:
        subprocess.run(['awk', '{print $2, $5, $6}', bfile1 + '_flipped.bim'],
                       stdout=out_file)
    
    # 3. remove uncorrected variants    
    with open(bfile1 + '_uncorrected_vars.tmp', 'w') as out_file:
        sort_proc = subprocess.Popen(['sort', bfile + '_vars.tmp', 
                                      bfile1 + '_flipped_vars.tmp'],
                                     stdout=subprocess.PIPE)
        
        uniq_proc = subprocess.Popen(['uniq', '-u'],
                                     stdin=sort_proc.stdout,
                                     stdout=subprocess.PIPE)
        
        awk_proc = subprocess.Popen(['awk', '{print $1}'],
                                    stdin=uniq_proc.stdout,
                                    stdout=subprocess.PIPE)
        
        sort_proc1 = subprocess.Popen(['sort', '-u'],
                                      stdin=awk_proc.stdout,
                                      stdout=out_file)
                                      
        sort_proc1.wait()
    
    subprocess.run(['plink', '--bfile', bfile, '--exclude', 
                    bfile1 + '_uncorrected_vars.tmp', '--make-bed', 
                    '--out', bfile + suffix])
    
    subprocess.run(['plink', '--bfile', bfile1 + '_flipped', '--exclude', 
                    bfile1 + '_uncorrected_vars.tmp', '--make-bed', 
                    '--out', bfile1 + suffix])

def merge(bfile, bfile1, out):
    sys.stdout = sys.__stdout__
    
    sfx = "_noA1A2.bim"
    strip_a1a2(bfile, sfx)
    strip_a1a2(bfile1, sfx)
    
    sfx1 = '_biallelic'
    remove_dup_rsids(bfile, sfx, sfx1)
    remove_dup_rsids(bfile1, sfx, sfx1)

    sfx2 = '_vars.tmp'
    print_variants(bfile, sfx1 + '.bim', sfx2)

    sfx3 = '_extract_vars'
    extract_variants(tmp_file(bfile1, sfx1),
                     tmp_file(bfile, sfx2),
                     sfx3)

    print_variants(bfile1, sfx1 + sfx3 + '.bim', sfx2)
    extract_variants(tmp_file(bfile, sfx1),
                     tmp_file(bfile1, sfx2),
                     sfx3)
    
    sfx4 = '_update_map'
    sync_maps(tmp_file(bfile, sfx1 + sfx3),
              tmp_file(bfile1, sfx1 + sfx3),
              sfx4)
    
    sfx5 = '_synced'
    sync_alleles(tmp_file(bfile, sfx1 + sfx3),
                 tmp_file(bfile1, sfx1 + sfx3 + sfx4),
                 sfx5)
    
    subprocess.run(['plink', '--bfile', tmp_file(bfile, sfx1 + sfx3 + sfx5),
                    '--bmerge', tmp_file(bfile1, sfx1 + sfx3 + sfx4 + sfx5),
                    '--make-bed', '--out', out])
    
    print("Done! Filesets are successfully merged.\n")

if __name__ == "__main__":
    
    import os
    import sys
    import subprocess
    import argparse
    
    cwd = os.getcwd()
    rel_tmp_path = "/tmp"
    tmp_path = cwd + rel_tmp_path
    
    if os.path.exists(tmp_path) == False:
        os.mkdir(tmp_path)
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--bfile", '-b', required=True, 
                        help = "plink binary file prefix")
    parser.add_argument("--bfile1", '-b1', required=True, 
                        help = "plink binary file prefix")
    parser.add_argument("--out", '-o', required=True,
                        help = "path and prefix for output fileset")
    
    args = parser.parse_args()
    merge(args.bfile, args.bfile1, args.out)
   
    
