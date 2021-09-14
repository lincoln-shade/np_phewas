#============================================
# Writes a vector of 95% CIs for odds ratios
#============================================

library(pacman)
p_load(data.table, magrittr)
make_or_95_ci <- function(or, l95, u95, or_ref, round_digits=2, flip_less_than_1=TRUE) {
  if (is.na(or)) {
    output <- paste0(" - ")
    return(output)
  } 
  if (flip_less_than_1) {
    output <- if (or_ref >= 1) {
      paste0(round(or, round_digits), " [", round(l95, round_digits), "-", round(u95, round_digits) , "]")
    } else if (or_ref < 1 & or_ref > 0) {
      paste0(round(1 / or, round_digits), " [", round(1 / u95, round_digits), "-", round(1 / l95, round_digits), "]")
    } else {
      paste0("Error")
    }
    return(output)
  } else {
    output <- ifelse(or > 0,
                     paste0(round(or, round_digits), " [", round(l95, round_digits), "-", round(u95, round_digits), "]"),
                     paste0("Error"))
    return(output)
  }
}
make_or_95_ci <- Vectorize(make_or_95_ci, vectorize.args = c("or", "l95", "u95", "or_ref"))
