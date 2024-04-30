#==========================================
# create linker file for rsids
#==========================================

library(data.table)
library(stringi)
bim <- fread("tmp/act_no_dup.bim")
rsid <- fread("raw_data/common_snps.txt", skip = 2)
rsid[, V1 := as.integer(stri_replace_first_fixed(V1, "chr", ""))]
rsid[, V2 := V2 + 1]
rsid[V5 == 1, V6 := stri_replace_first_fixed(V6, ",", "")]
alt_a1 <- stri_split_fixed(rsid[V5 == 2, V6], ",")
alt_a1 <- as.data.table(matrix(unlist(alt_a1), ncol = 3, byrow = TRUE))
identical(alt_a1[, paste0(V1, ",", V2, ",")], rsid[V5 == 2, V6])
rsid[V5 == 2, V6 := ..alt_a1$V1]
rsid[V5 == 2, V7 := ..alt_a1$V2]

rsid_linker <- merge(bim, 
                     rsid[V5 %in% 1:2], 
                     by.x = c("V1", "V4"), 
                     by.y = c("V1", "V2"), 
                     all.x = TRUE, 
                     suffixes = c("bim", "rsid"))

rsid_linker_use <- rsid_linker[((V6bim == V4rsid) & (V5bim == V6rsid)) |
                                 ((V6bim == V6rsid) & (V5bim == V4rsid)) |
                                 ((V6bim == V4rsid) & (V5bim == V7)) |
                                 ((V6bim == V7) & (V5bim == V4rsid)) |
                                 ((V6bim == V7) & (V5bim == V6rsid)) |
                                 ((V6bim == V6rsid) & (V5bim == V7))
                                 , ]
rsid_linker_use[, V3rsid := paste0(V3rsid, ":", V6bim, ":", V5bim)]
fwrite(rsid_linker_use[!(duplicated(V2)), .(V2, V3rsid)], 
       file = "tmp/rsid_linker.tmp", 
       col.names = FALSE, 
       quote = FALSE, 
       sep = " "
       )
