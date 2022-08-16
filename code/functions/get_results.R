
get_results <- function(files) {
  require(data.table)
  nfiles <- length(files)
  file_list <- vector(mode = "list", length = nfiles)
  for (i in 1:nfiles) {
    file_list[[i]] <- fread(files[i])
  }
  
  results <- rbindlist(file_list)
  return(results)
}
