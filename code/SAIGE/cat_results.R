
library(data.table)

get_results <- function(files) {
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  
  results <- rbindlist(file_list)
  return(results)
}


get_file_paths <- function(prefix, pattern) {
  paste0(prefix, list.files(prefix, pattern))
}

output_prefix <- "output/gwas/mega/saige/"
hs_files <- get_file_paths(output_prefix, "hs_saige_results_chr[0-9]*.txt$")
microinf_files <- get_file_paths(output_prefix, 
                                 "microinf_saige_results_chr[0-9]*.txt$")
grossinf_files <- get_file_paths(output_prefix, 
                                 "grossinf_saige_results_chr[0-9]*.txt$")

hs <- get_results(hs_files)
microinf <- get_results(microinf_files)
grossinf <- get_results(grossinf_files)

fwrite(hs, 
       file = "output/gwas/mega/saige/hs_saige_results.txt",
       quote = FALSE, 
       sep = " ")
fwrite(microinf, 
       file = "output/gwas/mega/saige/microinf_saige_results.txt",
       quote = FALSE, 
       sep = " ")
fwrite(grossinf, 
       file = "output/gwas/mega/saige/grossinf_saige_results.txt",
       quote = FALSE, 
       sep = " ")
