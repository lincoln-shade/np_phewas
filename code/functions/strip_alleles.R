#==========================================
# strips A1:A2 suffix from variant IDs
#==========================================


require(stringi)
strip_alleles <- function(x) {
  x <- stringi::stri_replace_last_regex(x, ":[ACTG]*:[ACTG]*", "")
}