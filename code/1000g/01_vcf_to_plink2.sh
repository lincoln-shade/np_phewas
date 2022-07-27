#!/bin/bash

#=========================================================================
# create plink filesets for each chromosome in 1000 genomes phase 3 data,
# then merge into one plink fileset
#=========================================================================

# create file to list fileset prefixes to merge
rm ./tmp/1000g_plink_fileset_prefixes.tmp

# create plink fileset for each chromosome.
# using --snps-only because some indels have multiple versions
# and plink can't account for this very well.
# have to write duplicate variant ids to file because some still present
# even with plink qc flags.

for chr in {1..22}
  do
  vcf=/data_global/1000g/hg38/ALL/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

  # # write duplicate variant ids to file
  # zcat $vcf | grep -v '^#' | cut -f 3 | sort | \
  #   uniq -d > tmp/1000g_chr"$chr"_dups.tmp

  plink2 \
    --vcf $vcf \
    --geno 0.01 \
    --make-pgen \
    --set-all-var-ids @:#_b38:\$r:\$a \
    --new-id-max-allele-len 51 \
    --vcf-half-call m \
    --maf 0.01 \
    --out tmp/1000g_chr"$chr".tmp

    echo tmp/1000g_chr"$chr".tmp >> tmp/1000g_plink_fileset_prefixes.tmp
  done

plink2 \
  --pmerge-list tmp/1000g_plink_fileset_prefixes.tmp \
  --make-pgen \
  --out data/1000g/1000g


# rm tmp/1000g*.tmp*

# plink2 \
#   --bfile data/1000g/1000g \
#   --make-pgen \
#   --out data/1000g/1000g