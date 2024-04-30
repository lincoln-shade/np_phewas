#============================================
# Writes a vector of 95% CIs for odds ratios
#============================================

make_or_95_ci <- function(or, l95, u95, or_ref, round_digits=2, 
                          flip_less_than_1=FALSE) {
  require(data.table)
  if (is.na(or)) {
    output <- paste0(" - ")
    return(output)
  } 
  if (flip_less_than_1) {
    output <- if (or_ref >= 1) {
      paste0(round(or, round_digits), " [", round(l95, round_digits), 
             "-", round(u95, round_digits) , "]")
    } else if (or_ref < 1 & or_ref > 0) {
      paste0(round(1 / or, round_digits), " [", round(1 / u95, round_digits), 
             "-", round(1 / l95, round_digits), "]")
    } else {
      paste0("Error")
    }
    return(output)
  } else {
    output <- ifelse(or > 0,
                     paste0(round(or, round_digits), " [", 
                            round(l95, round_digits), "-", 
                            round(u95, round_digits), "]"),
                     paste0("Error"))
    return(output)
  }
}
make_or_95_ci <- Vectorize(make_or_95_ci, 
                           vectorize.args = c("or", "l95", "u95", "or_ref"))

make_95CI <- function(or, p, ci=0.95) {
  beta <- log(or)
  lt_zero <- beta < 0
  z <- qnorm(p / 2, lower.tail = lt_zero)
  sd <- beta / z
  zrad <- qnorm((1 - ci) / 2, lower.tail = FALSE)
  upper <- round(exp(beta + zrad * sd), 2)
  lower <- round(exp(beta - zrad * sd), 2)
  conf_int <- paste0(lower, "-", upper)
  return(conf_int)
}

make_95CI <- Vectorize(make_95CI)

estimate_se_from_beta_and_p <- function(beta, p) {
  if (p <= 0 | p >= 1) {return(NA)}
  z <- qnorm(p / 2, lower.tail = FALSE)
  sd <- abs(beta) / z
  return(sd)
}

estimate_se_from_beta_and_p = Vectorize(estimate_se_from_beta_and_p)
