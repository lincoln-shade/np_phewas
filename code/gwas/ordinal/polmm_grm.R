library(data.table)
library(magrittr)
library(POLMM)

bfile <- "tmp/mega_np_prune"
n_parts_grm <- 14
chrs <- 1:22
tmp_dir <- "tmp/"

for (chr in chrs) {
  for (part in 1:n_parts_grm) {
    getSparseGRMParallel(chr, part, n_parts_grm, bfile,
                         tempDir = tmp_dir,
                         threadNum = 64)
  }
}

SparseGRM = getSparseGRM(chr, bfile, n_parts_grm, tempDir = tmp_dir)

bfile_name <- basename(bfile)
search_str <- paste0("Plink-", bfile_name, "-chr")
tmp_files <- list.files(tmp_dir) %>% 
  .[grep(search_str, .)] %>% 
  paste0("tmp/", .)

file.remove(tmp_files)

saveRDS(SparseGRM, file = "tmp/SparseGRM.Rds")
