# create file to list fileset prefixes to merge
echo > ./tmp/mc_caa_plink_fileset_prefixes.tmp

# create plink fileset for each chromosome.
# using --snps-only because some indels have multiple versions
# and plink can't account for this very well.
# have to write duplicate variant ids to file because some still present
# even with plink qc flags.

for batch in A B
  do
  for chr in 17 19
    do
    vcf=./data/mc_caa/Genetic_Variants/Imputed_Genotypes/MC-CAA_Batch"$batch"_chr"$chr".dose.vcf.gz

    # write duplicate variant ids to file
    # zcat $vcf | grep -v '^#' | cut -f 3 | sort | \
    # uniq -d > tmp/1000g_chr"$chr"_dups.tmp
    
    # plink \
    #   --vcf-filter \
    #   --vcf $vcf \
    #   --make-bed \
    #   --geno \
    #   --maf 0.01 \
    #   --biallelic-only strict \
    #   --snps-only just-acgt \
    #   --set-missing-var-ids @:#[b38]:\$1:\$2 \
    #   --vcf-half-call 'missing' \
    #   --extract-if-info R2 >= 0.8 \
    #   --out data/mc_caa/plink/mc_caa_batch"$batch"_chr"$chr"
    
    # echo tmp/1000g_chr"$chr".tmp >> tmp/1000g_plink_fileset_prefixes.tmp
    done
    
    plink \
      --bfile data/mc_caa/plink/mc_caa_batch"$batch"_chr17 \
      --recode A \
      --snp 17:8833591 \
      --out data/mc_caa/plink/mc_caa_batch"$batch"_chr17_rs72844606
    
    plink \
      --bfile data/mc_caa/plink/mc_caa_batch"$batch"_chr19 \
      --recode A \
      --snp 19:45454766 \
      --out data/mc_caa/plink/mc_caa_batch"$batch"_chr19_rs7247551
  done