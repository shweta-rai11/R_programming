library(affy)
library(GEOquery)
library(tidyverse)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("dplyr")
BiocManager::install("magrittr")
BiocManager::install("GEOquery")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("dplyr", "magrittr", "GEOquery", "vsn"))
BiocManager::install(c("DESeq2"))

library(affy)
library(GEOquery)
library(dplyr)
library(magrittr)
library(vsn)
# Set working directory
setwd("~/sweta/raw_files/Rhematoid_arthritis/UncompressedCEL")

# Reading files 
raw.data <- ReadAffy(celfile.path = "~/sweta/raw_files/Rhematoid_arthritis/UncompressedCEL")

#RMA normalization
normalized.data <- rma(raw.data)

# expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

# Map probe IDs to gene symbols
'''gse <- getGEO("GSE13837", GSEMatrix = TRUE)
feature.data <- gse$GSE148537_series_matrix.txt.gz@featureData@data
feature.data <- feature.data[, c(1, 11)]'''

# Load required library
library(GEOquery)

# Load the GEO dataset
gse <- getGEO("GSE13837")

# Extract the expression data
expression_data <- exprs(gse[[1]])

# Extract feature data
feature_data <- featureData(gse[[1]])

# Access the feature data
feature_data <- feature_data@data


normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')

# Function to clean data (handle NaN values and convert to numeric)
clean_data <- function(data) {
  cleaned_data <- apply(data, 2, function(x) ifelse(is.nan(x), median(x, na.rm = TRUE), x))
  cleaned_data_numeric <- apply(cleaned_data, 2, function(x) {
    x_numeric <- suppressWarnings(as.numeric(x))
    ifelse(is.na(x_numeric), median(x_numeric, na.rm = TRUE), x_numeric)
  })
  cleaned_data_df <- as.data.frame(cleaned_data_numeric)
  return(cleaned_data_df)
}

# Clean data
cleaned_normalized_expr <- clean_data(normalized.expr[, -1])
cleaned_normalized_data <- clean_data(exprs(normalized.data))

# Create boxplots
#par(mfrow = c(1, 2)) # Set up a 1x2 grid of plots
#boxplot(cleaned_normalized_expr, main = "Boxplot Before Normalization")
#boxplot(cleaned_normalized_data, main = "Boxplot After Normalization")


print(normalized.data)
str(normalized.data)
str(normalized.expr)
any(is.nan(exprs(normalized.data)))

# Histograms to compare before and after normalization with color
par(mfrow = c(1, 2))
hist(exprs(raw.data), col = "skyblue", main = "Histogram Before Normalization")
hist(exprs(normalized.data), col = "lightgreen", main = "Histogram After Normalization")



#plotting
# Open a PNG graphics device to save the plots
png("histograms.png")

# Set up the layout for side-by-side plots
par(mfrow = c(1, 2))

# Plot the histogram before normalization
hist(exprs(raw.data), col = "skyblue", main = "Histogram Before Normalization")

# Plot the histogram after normalization
hist(exprs(normalized.data), col = "lightgreen", main = "Histogram After Normalization")

# Close the PNG graphics device
dev.off()

