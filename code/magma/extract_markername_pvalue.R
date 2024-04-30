#!/usr/bin/env Rscript

# Load necessary libraries
library(argparse)
library(data.table)

# Setting up Argument Parser
parser <- ArgumentParser()

# Adding arguments
parser$add_argument("-i", "--input", type = "character", help = "Input CSV file name")
parser$add_argument("-s", "--snp_column", type = "character", default = "MarkerName", help = "Column name for SNP, default is 'MarkerName'")
parser$add_argument("-p", "--pvalue_column", type = "character", default = "P.value", help = "Column name for P-value, default is 'P.value'")
parser$add_argument("-o", "--output", type = "character", help = "Output CSV file name")

# Parse the arguments
args <- parser$parse_args()

# Read the input file
data <- fread(args$input)

# Check if the SNP and P-value columns exist
if (!(args$snp_column %in% colnames(data)) || !(args$pvalue_column %in% colnames(data))) {
    stop("Specified columns not found in the data")
}

# Extract SNP and P-value columns
output_data <- data[, .(SNP=get(..args$snp_column), P=get(..args$pvalue_column))]

# Write the extracted data to a new CSV file
fwrite(output_data, file = args$output, sep = " ")

cat("File saved:", args$output, "\n")
