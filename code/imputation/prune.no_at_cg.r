options(stringsAsFactors = FALSE)

list.of.packages <- c("dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(dplyr)

argv <- commandArgs(TRUE)


all.1000genomes = read.table(paste(argv[1], "/ALL.1000genemes.no-at-cg.uniq.rs", sep = ""), h = F)
# plink.prune.in = read.table(paste("tmp/", argv[2], ".prune.in", sep = ""), h = F)
plink.prune.in = read.table(paste("tmp/plink.prune.in", sep = ""), h = F)

int.rs = plink.prune.in %>% filter(V1 %in% all.1000genomes$V1)

write.table(int.rs, "tmp/prune.1000genomes.no_at_cg.in", quote = F, col.names = F, row.names = F)


