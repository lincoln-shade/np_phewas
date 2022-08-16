
#====================================================
# Get SAIGE results files by chromosom and rbindlist
#====================================================

library(data.table)
library(argparse)
library(stringi)

parser <- ArgumentParser()
parser$add_argument("-p", "--prefix", 
                    help = paste0("prefix for results to rbind. ",
                                  "Files should end in '_chr[0-9]*.txt'. ",
                                  "prefix should end in '_chr'."))
parser$add_argument("-o", "--output", 
                    help = paste0(
                      "path for output file. If not used, will ",
                      "default to --prefix minus '_chr' plus '.txt'"))
args <- parser$parse_args()
print(args)

get_results <- function(files) {
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  
  results <- rbindlist(file_list)
  return(results)
}

files <- list.files(dirname(args$prefix), 
                    paste0(basename(args$prefix), "[0-9]*.txt$"))
print(files)

results <- get_results(paste0(dirname(args$prefix), "/", files))
setorder(results, CHR, POS)

if (is.null(args$output)) {
  out <- stri_replace_last_fixed(args$prefix, "_chr", ".txt")
} else {
  out <- args$output
}

fwrite(results, file = out, sep = " ", quote = FALSE, na = "NA")