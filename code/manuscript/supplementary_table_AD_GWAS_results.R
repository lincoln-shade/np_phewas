
library(data.table)
ad_gwas <- 
  fread("/data_global/summary_statistics/GCST90027158_buildGRCh38.tsv.gz")

top_snps <- c("rs6733839", "rs6460900", "rs4721058", "rs2000660", "rs72807981",
              "rs769449", "rs12721051", "rs4420638", "rs4803778"
)

ad_gwas_top_snps <- ad_gwas[variant_id %in% top_snps]

jansen_2019 <- fread("/data_global/summary_statistics/jansenIE2019_AD/AD_sumstats_Jansenetal_2019sept.txt.gz")

jansen_2019_top_snps <- jansen_2019[SNP %in% top_snps]
jansen_2019_top_snps[, chr_bp := paste0(CHR, "_", BP)]

wightman <- fread("/data_global/summary_statistics/wightman2021_AD/PGCALZ2sumstatsExcluding23andMe.txt.gz")
wightman[, CHR_BP := paste0(chr, "_", PosGRCh37)]
wightman_top_snps <- wightman[CHR_BP %in% jansen_2019_top_snps$chr_bp]
