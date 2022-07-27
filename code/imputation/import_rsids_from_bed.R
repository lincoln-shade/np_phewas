#==========================================
# create linker file for rsids
#==========================================

library(data.table)
library(stringi)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--pvar", help = "PLINK2 .pvar file")
parser$add_argument("--bed", help = "BED file with variant info",
                    default = "raw_data/dbSnp153Common.bed")
parser$add_argument("--out", help = "path to output linker file")
args <- parser$parse_args()

pvar <- fread(args$pvar)
rsid <- fread(cmd = paste("awk '{print $1,$2,$3,$4,$5,$6,$7}'", args$bed))
rsid[, V1 := as.integer(stri_replace_first_fixed(V1, "chr", ""))]
rsid[, V2 := V2 + 1L]
rsid[V6 == "1", V7 := stri_replace_first_fixed(V7, ",", "")]
alt_a1 <- stri_split_fixed(rsid[V6 == "2", V7], ",")
alt_a1 <- as.data.table(matrix(unlist(alt_a1), ncol = 3, byrow = TRUE))
identical(alt_a1[, paste0(V1, ",", V2, ",")], rsid[V6 == 2, V7])
rsid[V6 == "2", V7 := ..alt_a1$V1]
rsid[V6 == "2", V8 := ..alt_a1$V2]
rsid <- rsid[V6 %in% c("1", "2")]

# merge by Chromosome and base pair
rsid_linker <- merge(pvar, 
                     rsid, 
                     by.x = c("#CHROM", "POS"), 
                     by.y = c("V1", "V3"), 
                     all.x = TRUE, 
                     suffixes = c("pvar", "rsid"))

# keep only variants where alleles match
# (e.g. REF == A1 and ALT == A2 or REF == A2 and ALT == A1)
rsid_linker_use <- rsid_linker[((REF == V5) & (ALT == V7)) |
                                 ((REF == V7) & (ALT == V5)) |
                                 ((REF == V5) & (ALT == V8)) |
                                 ((REF == V8) & (ALT == V5))]
rsid_linker_use[V4 %in% V4[duplicated(V4)], 
                V4 := paste0(V4, ":", REF, ":", ALT)]
fwrite(rsid_linker_use[!(duplicated(ID)), .(ID, V4)], 
       file = args$out, 
       col.names = FALSE, 
       quote = FALSE, 
       sep = " "
)
