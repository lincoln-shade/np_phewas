# format METAL results for downstream work

#!/usr/bin/env Rscript --vanilla

# Load the required packages
library(data.table)
library(stringi)
library(argparse)

# function for flipping direction of effect in Direction column (e.g. "++-" -> "--+")
flip_effect_directions_fast <- function(directions_vector) {
  # Temporarily replace "+" with a placeholder ("p")
  temp_vector <- stri_replace_all_regex(directions_vector, "\\+", "p")
  
  # Replace "-" with "+"
  temp_vector <- stri_replace_all_regex(temp_vector, "-", "\\+")
  
  # Replace placeholder ("p") with "-"
  flipped_vector <- stri_replace_all_regex(temp_vector, "p", "-")
  
  return(flipped_vector)
}

# function for flipping alleles, direction of effect, allele frequency, and direction for data table
flip_alleles_effect_direction <- function(dt) {
  setnames(merged_dt, c("Allele1", "Allele2"), c("Allele2", "Allele1"))
  dt[, Freq1 := 1 - Freq1]
  dt[, Effect := -Effect]
  dt[, MinFreq1 := 1 - MaxFreq]
  dt[, MaxFreq := 1 - MinFreq]
  dt[, MinFreq := NULL]
  setnames(dt, "MinFreq1", "MinFreq")
  dt[, Direction := flip_effect_directions_fast(Direction)]
  return(dt)
}
# Define parser for command-line arguments
parser <- ArgumentParser(description = "Merge METAL meta-analysis GWAS file with PLINK .bim file")

# Add arguments to the parser
parser$add_argument("--metal_file", required = TRUE, help = "Path to METAL meta-analysis GWAS file")
parser$add_argument("--bim_file", required = TRUE, help = "Path to PLINK .bim file")
parser$add_argument("--output_file", required = TRUE, help = "Path to output merged file")
parser$add_argument("--flip_alleles", default = "FALSE", help = "Flip alleles and effect sizes? Default = FALSE")

# Parse the command-line arguments
args <- parser$parse_args()
args$flip_alleles = as.logical(args$flip_alleles)

# Read METAL meta-analysis GWAS file into a data.table
# Assume columns: MarkerName, P-value (add more as needed)
metal_dt <- fread(args$metal_file, check.names = TRUE)

# Read PLINK .bim file into a data.table
# Columns: Chromosome, MarkerName, GeneticDistance, BasepairPosition
bim_dt <- fread(
  args$bim_file, 
  col.names = c("chr", "marker", "gd", "bp", "a1", "a2")
)

# Merge the two data.tables based on the MarkerName
merged_dt <- merge(
  metal_dt, 
  bim_dt[, .(marker, chr, bp)], 
  by.x = "MarkerName", 
  by.y = "marker", 
  all.x = TRUE
)

# capitalize the allele values
merged_dt[, Allele1 := toupper(Allele1)]
merged_dt[, Allele2 := toupper(Allele2)]

# Make Odds ratios and allele frequencies with respect to the minor allele
merged_dt_freq_lt50 = merged_dt[Freq1 <= 0.5]
merged_dt_freq_gt50 = flip_alleles_effect_direction(merged_dt[Freq1 > 0.50])
rm(merged_dt)
merged_dt1 = rbindlist(list(merged_dt_freq_lt50, merged_dt_freq_gt50), use.names = TRUE)


# set columns in defined order and order by chromosome and basepair
setcolorder(merged_dt1, c("chr", "bp", "MarkerName", "Allele1", "Allele2", "Freq1",	"FreqSE",	"MinFreq",	"MaxFreq",	"Effect",	"StdErr",	"P.value", "Direction"))
setorder(merged_dt1, chr, bp)
merged_dt1 = merged_dt1[MinFreq >= 0.01] # only keep variants with MAF >= 1% in each data source used

# Write the merged data.table to the specified output file
fwrite(merged_dt1, args$output_file)

# Print out a message indicating that the merge was successful
cat("Successfully merged", args$metal_file, "and", args$bim_file, "into", args$output_file, "\n")

