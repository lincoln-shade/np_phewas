# UKY ADC-imputation
Path to files: /data_global/UKY_ADC/TOPMedImputation/ \
TOPMed Imputation Server: https://imputation.biodatacatalyst.nhlbi.nih.gov/#!

### hg19
#### 0. Preparation (run just once) 
```
(1) Download ALL.TOPMed_freeze3_hg19_dbSNP.vcf.gz from https://bravo.sph.umich.edu/freeze3a/hg19/download
(2) Run ./CreateTOPMed.pl -i ALL.TOPMed_freeze3_hg19_dbSNP.vcf.gz
(3) Download perl script from https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip 
```

#### 1. Do QC
```
$ python qc.py --bfile PLINK_FILE_NAME

Options: 
 --mind 0.05 (missing rate to remove individuals) 
 --maf 0.05 (MAF to remove rare variants)
 --geno 0.05 (missing rate to remove variants)
 --ibd 0.185 (PHI_HAT to remove one of pairs)
 --keep FILE (filename with FID and IID to extract a subset of individuals)
 --pheno FILE (phenotype file with FID, IID, binary phenotype 0 = unaffected, 1 = affected)
 --pheno-name PHENOTYPE_COLUMN_NAME (phenotype name in phenotype file)
 --out FILE (filename for output)

Output: 
  PLINK_FILE_NAME.qced.bim/fam/bed (QCed PLINK filesets) 

Path to qc.py:
  /data_global/code/qc/qc.py
```
``` diff
- PTID=451 and 778 are in both adgcr519 and adgcr9_19
```
#### 2. Generate a file with FID and IID to keep or remove a subset of individuals using PCs and PCA plot created in 1.

#### 3. Prepare .vcf.gz files in each chromosome for TOPMed imputation 
```
$ python preparation.py --bfile PLINK_FILE_NAME.qced

Options:
 --keep FILE (filename with FID and IID to keep a subset of individuals)
 --remove FILE (filename with FID and IID to remove a subset of individuals)

Output: 
  PLINK_FILE_NAME_hwe_extract-updated-chr1.vcf.gz 
  .... 
  PLINK_FILE_NAME_hwe_extract-updated-chr22.vcf.gz 
```

#### 4. Upload .vcf.gz files on TOPMed Imputation Server
```
https://imputation.biodatacatalyst.nhlbi.nih.gov/#!
```

#### 5. Set parameters on TOPMed Imputation Server
![alt text](https://github.com/kyka222/ADC-imputation/blob/master/TopMedImputation.setting.PNG?raw=true=50x50)

#### 6. Download imputed chr_1.zip to chr_22.zip on statgen
```diff
- zipped vcf files created with GRCh38 (hg38)
```

#### 7. Unzip with password (which is sent by email) in adc1 to adc10
```
$ seq 1 22 | xargs -n 1 bash -c 'unzip -P "PASSWORD" chr_$0.zip'
```
#### 8. Create index for each .dose.vcf.gz in adc1 to adc10
```
$ seq 1 22 | xargs -n 1 bash -c 'tabix -p vcf chr$0.dose.vcf.gz'
```
#### 9. Merge adc1 to adc10 for each chr (Error: IID=451 and 778 are duplicated in adc5 and adc9 --> add --force-samples)
```
$ seq 1 22 | xargs -n 1 bash -c 'bcftools merge adc*/chr$0.dose.vcf.gz --force-samples -Oz -o merge/chr$0.alladcs.dose.vcf.gz' 
```
#### 10. Convert vcf to plink
``` diff
$ seq 1 22 | xargs -n 1 bash -c 'plink --vcf ../vcf/chr$0.alladcs.dose.vcf.gz --make-bed --out chr$0.alladcs'

- preQC plink files could not be merged due to memory limitation. 
- Therefore, varinats with maf < 0.01 and/or missing rate > 0.05 were removed before merge
- rs ID looks like chr1:50001586:G:A ---> need to convert it to rs ID based on chromosome and position
```
***
#### 11. remove variants with maf < 0.01 and/or missing rate > 0.05 
```
$ seq 1 22 | xargs -n 1 -P 22 bash -c 'plink --bfile preQC/chr$0.alladcs --maf 0.01 --geno 0.05 --make-bed --out postQC/chr$0.alladcs'
```
#### 12. merge chr1 to chr22
```
$ seq 2 22 | xargs -n 1 bash -c 'echo chr$0.alladcs >> merge.list'
$ plink --bfile chr1.alladcs --merge-list merge.list --make-bed --out alladcs
```
#### 13. QC
```
$ python qc.py --bfile alladcs --maf 0.01 --geno 0.05 --out alladcs.qced
```

***
#### 14. Convert chr1:50001586:G:A to rs ID
(1) access UCSC mySQL https://genome.ucsc.edu/goldenPath/help/mysql.html
```
$ mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -P 3306
```
(2) get database and table names
```
mysql> show databases;
mysql> use hg38;
mysql> show tables;
mysql> exit
```
(3) get rs IDs
```
$ seq 1 22 | xargs -n 1 -P 5 bash -c 'mysql -h genome-mysql.soe.ucsc.edu -u genome -A -P 3306 -D hg38 --skip-column-names -e "select * from snp151 where chrom = \"chr$0\";" > ucsc.hg38.snp151.chr$0.rs'
```
(4) create linker file with chr position rs major minor
```
$ seq 1 22 | xargs -n 1 -P 5 bash -c 'cat ucsc.hg38.snp151.chr$0.rs | awk '\''{print $2,$3+1,$5,$10}'\'' | grep -v lengthTooLong | sed -e '\''s/chr//g'\'' | sed -e '\''s/\// /g'\'' > ucsc.hg38.snp151.chr$0.linker'
```
(5) create files with old rs and new rs for PLINK --update-name option (see http://www.cog-genomics.org/plink/1.9/data)
```
$ seq 1 22 | xargs -n 1 -P 22 bash -c 'sh convert.hg38.rs.sh $0' ## DON'T RUN! It takes long time. All files are /data_global/UKY_ADC/TOPMedImputation/linker
```
(6) combine all linkers
```
$ seq 1 22  | xargs -n 1 bash -c 'cat chr$0.rs.linker >> chrall.rs.linker'
```
(7) update rs ID
```
$ plink --bfile alladcs.qced --update-name chrall.rs.linker --make-bed --out alladcs.qced.rs
```

***
#### 15. Run pca
```
$ python pca.py --bfile PLINK_FILE_NAME

Options: 
 --superpop ALL (choose superpopulation. As of 7/31/2020, only ALL is available) 
 --threads 1 (# of threads used to convert 1000g vcf to plink files)
 --build hg19 (hg18 or hg19 or hg38)
 --out FILE (filename for output)
 --prune-in FILE_NAME (prune.in file if you have. If not specified, prune.in file will be created with your data)
 --umap (whether umap is performed. Default = False)

Output: 
  PLINK_FILE_NAME.pcs-ancestry-plot.pdf
  PLINK_FILE_NAME.pcs.csv (PC1 to PC20)

Path to pca.py:
  /data_global/code/pca/pca.py
  ```
