#! /usr/bin/Rscript --vanilla

library(data.table)
library(snpStats)
library(stringi)
library(foreach)
library(doParallel)
library(argparse)
library(stringi)

parser <- ArgumentParser()
parser$add_argument("-b", "--bim", help = "PLINK .bim file",
                    default = "data/mega/mega_np.bim")
parser$add_argument("-m", "--snp_set_size", 
                    help = "max number of SNPs per SnpMatrix",
                    default = 50000L)
parser$add_argument("-o", "--out", default = "tmp/", 
                    help = "prefix for output SnpMatrix .Rds objects")
parser$add_argument("-c", "--cores", default = 64, 
                    help = "number of cores for parallel computing")
args <- parser$parse_args()

if (is.null(args$bim) | is.null(args$snp_set_size)) {
  stop("BIM file and SNP set size are required arguments.")
}

if (!dir.exists(args$out)) {
  dir.create(args$out)
}

plink <- stri_replace_last_fixed(args$bim, ".bim", "")
m <- as.integer(args$snp_set_size)
cores <- as.integer(args$cores)

# get total number of SNPs in .bim file
n_snp <- as.integer(
  stri_split_fixed(
    system(paste0("wc -l ", plink, ".bim"), intern = TRUE), 
    pattern=" ")[[1]][1]
)

n_groups <- (n_snp %/% m)
groups <- c(seq(0, m * n_groups, m), n_snp)

# use multicore, set to the number of our cores
registerDoParallel(cores = cores)  
out_list <- 
  foreach (i=2:length(groups), .combine = c) %dopar% {
    snp_start <- groups[i - 1] + 1
    snp_stop <- groups[i]
    mat <- read.plink(bed = paste0(plink, '.bed'), 
                      bim = paste0(plink, '.bim'), 
                      fam = paste0(plink, '.fam'),
                      select.snps = snp_start:snp_stop)
    out_file <- paste0(args$out, "SnpMatrix_", i - 1,".Rds")
    saveRDS(mat, file = out_file)
    rm(mat)
    out_file
  }
stopImplicitCluster()

saveRDS(out_list, file = paste0(args$out, "SnpMatrix_list.Rds"))


