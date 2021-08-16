#!/bin/bash

#=========================================================================
# create plink filesets for each chromosome in 1000 genomes phase 3 data,
# then merge into one plink fileset
#=========================================================================

# create file to list fileset prefixes to merge
echo > ./data/tmp/1000g_plink_fileset_prefixes.tmp

# create plink fileset for each chromosome.
# using --snps-only because some indels have multiple versions
# and plink can't account for this very well.
# have to write duplicate variant ids to file because some still present
# even with plink qc flags.

for chr in {1..22}
  do
  vcf=/data_global/1000g/hg38/ALL/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

  # write duplicate variant ids to file
  zcat $vcf | grep -v '^#' | cut -f 3 | sort | uniq -d > ./data/tmp/1000g_chr"$chr"_dups.tmp

  plink \
    --vcf-filter \
    --vcf $vcf \
    --make-bed \
    --geno \
    --maf 0.01 \
    --biallelic-only strict \
    --snps-only just-acgt \
    --set-missing-var-ids @:#[b37]:\$1:\$2 \
    --vcf-half-call 'missing' \
    --exclude ./data/tmp/1000g_chr"$chr"_dups.tmp \
    --out ./data/tmp/1000g_chr"$chr".tmp

    echo ./data/tmp/1000g_chr"$chr".tmp >> ./data/tmp/1000g_plink_fileset_prefixes.tmp
  done
  
plink \
  --merge-list ./data/tmp/1000g_plink_fileset_prefixes.tmp \
  --make-bed \
  --out data/1000g/1000g

rm ./data/tmp/1000g*.tmp*