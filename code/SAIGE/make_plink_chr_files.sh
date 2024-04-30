
# mkdir data/adc/plink
# # make chromosome-specific plink filesets
# for chr in $(seq 1 22)
#   do 
#   plink \
#     --bfile data/adc/adc_np \
#     --chr "$chr" \
#     --make-bed \
#     --out data/adc/plink/adc_np_chr"$chr"
#   done

# mkdir data/rosmap/plink
# # make chromosome-specific plink filesets
# for chr in $(seq 1 22)
#   do 
#   plink \
#     --bfile data/rosmap/rosmap_np \
#     --chr "$chr" \
#     --make-bed \
#     --out data/rosmap/plink/rosmap_np_chr"$chr"
#   done

mkdir data/act/plink
# make chromosome-specific plink filesets
for chr in $(seq 1 22)
  do 
  plink \
    --bfile data/act/act_np \
    --chr "$chr" \
    --make-bed \
    --out data/act/plink/act_np_chr"$chr"
  done