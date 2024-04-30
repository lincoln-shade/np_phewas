#!/usr/bin/env Rscript

# Load the data.table package
library(data.table)

# Function to count non-missing rows in a specified column
count_non_missing <- function(file_name, column_name) {
    # Read the file
    data <- fread(file_name)
    
    # Check if the column exists
    if (!column_name %in% colnames(data)) {
        stop("Column not found in the data")
    }
    
    # Count non-missing (non-NA) values in the specified column
    count <- sum(!is.na(data[[column_name]]))
    
    return(count)
}

# Main script
args <- commandArgs(trailingOnly = TRUE)

# Check if there are at least two arguments
if (length(args) < 2) {
    stop("Two arguments required: file name and column name")
}

# Assign arguments to variables
file_name <- args[1]
column_name <- args[2]

# Call the function and print the result
result <- count_non_missing(file_name, column_name)
cat(result)