
library(data.table)
get_results_list <- function(files) {
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  return(file_list)
}