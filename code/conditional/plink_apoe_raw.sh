
#======================================================
# get APOE haplotype SNP info for participants.
# APOE info at https://www.snpedia.com/index.php/APOE
#======================================================

rs429358_pos=44908684
rs7412_pos=44908822

cat data/adc/adc.bim | \
  awk -v bp1=$rs429358_pos -v bp2=$rs7412_pos \
    '$1 == 19 && ($4 == bp1 || $4 == bp2) {print $2}'

find_apoe_varids () {
  cat $1 | \
  awk -v bp1=$rs429358_pos -v bp2=$rs7412_pos \
    '$1 == 19 && ($4 == bp1 || $4 == bp2) {print $2}'
}

# NACC
find_apoe_varids data/adc/adc.bim > tmp/adc_apoe_varids.tmp

plink \
  --bfile data/adc/adc \
  --chr 19 \
  --extract tmp/adc_apoe_varids.tmp \
  --recode A \
  --out data/adc/apoe_snps

# ROSMAP
rosmap_prefix=/data_global/ROSMAP/rosmap_topmed_20210303/ROSMAP_NHW_imputed_final
find_apoe_varids "$rosmap_prefix".bim > tmp/rosmap_apoe_varids.tmp

plink \
  --bfile "$rosmap_prefix" \
  --fam "$rosmap_prefix"_converted.fam \
  --chr 19 \
  --extract tmp/rosmap_apoe_varids.tmp \
  --recode A \
  --out data/rosmap/apoe_snps

# ACT
act_prefix=data/act/act
find_apoe_varids "$act_prefix".bim > tmp/act_apoe_varids.tmp

plink \
  --bfile "$act_prefix" \
  --chr 19 \
  --extract tmp/act_apoe_varids.tmp \
  --recode A \
  --out data/act/apoe_snps
