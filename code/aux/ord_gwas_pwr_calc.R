#! /usr/bin/Rscript --vanilla

library(rms)
library(ggplot2)
library(data.table)
library(scales)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--odds_ratios", default = "c(1.35, 1.45, 1.55, 1.65)",
                    help = "string of vector of ORs to be parsed")
parser$add_argument("--sims", default = "5e4", 
                    help = "number of simulations to run")
parser$add_argument("--n", default = "3005", 
                    help = "number of subjects")
parser$add_argument("--alpha", default = "5e-8", 
                    help = "number of simulations to run")
parser$add_argument("--mafs", default = "seq(0.05, 0.5, 0.01)", 
                    help = "string of vector of MAFs to be parsed")
parser$add_argument("--out", help = ".Rds file path to save results")
args <- parser$parse_args()

odds_ratios <- eval(parse(text = args$odds_ratios))
sims <- as.numeric(args$sims)
n <- as.numeric(args$n)
alpha <- eval(parse(text = args$alpha))
mafs <- eval(parse(text = args$mafs))
sys_date <- Sys.Date()[1]

if (is.null(args$out)) {
  out <- paste0("ord_power_n", n, "_p", alpha, "_", sys_date, ".Rds")
} else if (grep(".Rds|.rds|.RDS", args$out) == 0) {
  stop("--out file needs to end in .Rds, .rds, or .RDS")
} else {
  out <- args$out
}

tmpfun <- function(n=2800, maf=0.5, or=1, beta0=-1, beta1=0.5, beta2=0.5) {
  require(rms)
  x <- rbinom(n, 2, maf)
  eta1 <- beta0 + log(or) * x
  eta2 <- eta1 + beta1
  eta3 <- eta2 + beta2
  p1 <- exp(eta1)/(1+exp(eta1))
  p2 <- exp(eta2)/(1+exp(eta2))
  p3 <- exp(eta3)/(1+exp(eta3))
  tmp <- runif(n)
  y <- (tmp < p1) + (tmp < p2) + (tmp < p3)
  table(y)
  fit <- lrm(y~x)
  fit$stats[5]
}

calc_power <- function(n, ors, mafs, sims, alpha) {
  dt <- data.table(MAF = mafs)
  for (i in 1:length(ors)) {
    power <- numeric(length(mafs))
    or_colname <- paste0(ors[i])
    for (j in 1:length(mafs)) {
      out <- replicate(sims, tmpfun(maf = mafs[j], or = ors[i]))
      power[j] <- mean( out < alpha )
    }
    dt[, (or_colname) := power[order(power)]]
    dtl <- melt(dt[], id.vars = 'MAF', 
                variable.name = 'OR', value.name = 'Power')
  }
  return(dtl)
}

dtl <- calc_power(ors = odds_ratios,
                  sims = sims,
                  mafs = mafs,
                  alpha = alpha)
saveRDS(dtl, file = out)
