
#============================
# merge two plink filesets
#============================

def tmp_file(pfile, suffix=''):
    return tmp_path + "/" + os.path.basename(pfile) + suffix

def strip_a1a2(pfile, suffix):
    # strips the A1:A2 suffixes from variant rsIDs
    out_path = tmp_file(pfile, suffix)
    with open(out_path, "w") as file_out:
        cat_proc = subprocess.Popen(["cat", pfile + ".pvar"], 
                                    stdout=subprocess.PIPE)
        awk_proc = subprocess.Popen(["awk", 
                                     "{print gensub(/:[ACTG]*:[ACTG]*/, \"\", 1)}"],
                                    stdin=cat_proc.stdout,
                                    stdout=file_out)
        awk_proc.wait()

def remove_dup_rsids(pfile, suffix, suffix1):
    bim_file = tmp_file(pfile, suffix)
    out_path = tmp_file(pfile, "_dup_rsids.tmp")
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
    
    cmd = ['plink', '--pfile', pfile, '--bim', bim_file, '--make-bed',
           '--out', tmp_file(pfile, suffix1)]

    if os.stat(out_path).st_size != 0:
        cmd += ['--exclude', out_path]

    subprocess.run(cmd)

def print_variants(pfile, suffix, suffix1):
    with open(tmp_file(pfile, suffix)) as in_file:
        with open(tmp_file(pfile, suffix1), 'w') as out_file:
            subprocess.run(['awk', '{print$2}'],
                           stdin=in_file,
                           stdout=out_file)

def extract_variants(pfile, variants, suffix):
    cmd = ['plink', '--pfile', pfile, '--extract', variants, 
           '--make-bed', '--out', tmp_file(pfile, suffix)]
    subprocess.run(cmd)
    
def sync_maps(pfile, pfile1, suffix):
    subprocess.run(['plink', '--pfile', pfile1, '--update-map', 
                    pfile + '.pvar', '4', '2', '--make-bed', 
                    '--out', pfile1 + suffix])

def sync_alleles(pfile, pfile1, suffix):
    # 1. set pfile1 alleles to pfile alleles
    subprocess.run(['plink', '--pfile', pfile1, '--a1-allele', pfile + '.pvar',
                    '5', '2', '--make-bed', '--out', pfile1 + '_upmap'])
    
    # 2. try to resolve strand issues
    with open(pfile + '_vars.tmp', 'w') as out_file:
        subprocess.run(['awk', '{print $2, $5, $6}', pfile + '.pvar'],
                       stdout=out_file)
    
    with open(pfile1 + '_vars.tmp', 'w') as out_file:
        subprocess.run(['awk', '{print $2, $5, $6}', pfile1 + '_upmap.pvar'],
                       stdout=out_file)
    
    with open(pfile1 + '_flip_list.tmp', 'w') as out_file: 
        sort_proc = subprocess.Popen(['sort', pfile + '_vars.tmp', 
                                      pfile1 + '_vars.tmp'], 
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
    
    subprocess.run(['plink', '--pfile', pfile1 + '_upmap', 
                    '--flip', pfile1 + '_flip_list.tmp',
                    '--a1-allele', pfile + '.pvar', '5', '2', '--make-bed',
                    '--out', pfile1 + '_flipped'])
    
    with open(pfile1 + '_flipped_vars.tmp', 'w') as out_file:
        subprocess.run(['awk', '{print $2, $5, $6}', pfile1 + '_flipped.pvar'],
                       stdout=out_file)
    
    # 3. remove uncorrected variants    
    with open(pfile1 + '_uncorrected_vars.tmp', 'w') as out_file:
        sort_proc = subprocess.Popen(['sort', pfile + '_vars.tmp', 
                                      pfile1 + '_flipped_vars.tmp'],
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
    
    subprocess.run(['plink', '--pfile', pfile, '--exclude', 
                    pfile1 + '_uncorrected_vars.tmp', '--make-bed', 
                    '--out', pfile + suffix])
    
    subprocess.run(['plink', '--pfile', pfile1 + '_flipped', '--exclude', 
                    pfile1 + '_uncorrected_vars.tmp', '--make-bed', 
                    '--out', pfile1 + suffix])

def merge(pfile, pfile1, out):
    sys.stdout = sys.__stdout__
    
    sfx = "_noA1A2.pvar"
    strip_a1a2(pfile, sfx)
    strip_a1a2(pfile1, sfx)
    
    sfx1 = '_biallelic'
    remove_dup_rsids(pfile, sfx, sfx1)
    remove_dup_rsids(pfile1, sfx, sfx1)

    sfx2 = '_vars.tmp'
    print_variants(pfile, sfx1 + '.pvar', sfx2)

    sfx3 = '_extract_vars'
    extract_variants(tmp_file(pfile1, sfx1),
                     tmp_file(pfile, sfx2),
                     sfx3)

    print_variants(pfile1, sfx1 + sfx3 + '.pvar', sfx2)
    extract_variants(tmp_file(pfile, sfx1),
                     tmp_file(pfile1, sfx2),
                     sfx3)
    
    sfx4 = '_update_map'
    sync_maps(tmp_file(pfile, sfx1 + sfx3),
              tmp_file(pfile1, sfx1 + sfx3),
              sfx4)
    
    sfx5 = '_synced'
    sync_alleles(tmp_file(pfile, sfx1 + sfx3),
                 tmp_file(pfile1, sfx1 + sfx3 + sfx4),
                 sfx5)
    
    subprocess.run(['plink', '--pfile', tmp_file(pfile, sfx1 + sfx3 + sfx5),
                    '--bmerge', tmp_file(pfile1, sfx1 + sfx3 + sfx4 + sfx5),
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
    
    parser.add_argument("--pfile", '-p', required=True, 
                        help = "plink binary file prefix")
    parser.add_argument("--pfile1", '-p1', required=True, 
                        help = "plink binary file prefix")
    parser.add_argument("--out", '-o', required=True,
                        help = "path and prefix for output fileset")
    
    args = parser.parse_args()
    merge(args.pfile, args.pfile1, args.out)
   
    
