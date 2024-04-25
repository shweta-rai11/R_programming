# Loading necessary libraries
library(affy)
library(GEOquery)
library(dplyr)
library(tidyr)  # For data manipulation functions
library(tibble)  # For rownames_to_column function

# Set working directory
setwd("~/sweta/raw_files/Rhematoid_arthritis")
# Set your working directory


# Create a new directory for .CEL files
dir.create("SeparatedCEL", showWarnings = FALSE)

# List all files in the current directory
files <- list.files()

# Filter out .CEL.gz files
cel_files <- files[grepl("\\.CEL\\.gz$", files)]

# Move .CEL.gz files to the SeparatedCEL directory
for (file in cel_files) {
  dest_file <- file.path("SeparatedCEL", basename(file))
  con <- gzfile(file, "rb") # Open .gz file
  out_con <- gzfile(dest_file, "wb") # Open destination .gz file
  while (length(data <- readLines(con, n = 1)) > 0) { # Reading line by line
    writeLines(data, con = out_con) # Write line to destination .gz file
  }
  close(con)
  close(out_con)
}

# Creating a new directory for uncompressed .CEL files
dir.create("UncompressedCEL", showWarnings = FALSE)

# Listing all files in the current directory
files <- list.files()

# Filtering out .CEL.gz files
cel_gz_files <- files[grepl("\\.CEL\\.gz$", files)]

# Unziping .CEL.gz files and save to the UncompressedCEL directory
for (file in cel_gz_files) {
  dest_file <- file.path("UncompressedCEL", sub("\\.gz$", "", basename(file))) # remove .gz extension from file name
  con <- gzfile(file, "rb") # Open .gz file
  out_con <- file(dest_file, "wb") # Open destination file without compression
  while (length(data <- readBin(con, "raw", n = 65536)) > 0) { # Read binary data in chunks
    writeBin(data, out_con) # Write binary data to destination file
  }
  close(con)
  close(out_con)
}


# Reading CEL files
raw.data <- ReadAffy(celfile.path = "~/sweta/raw_files/Rhematoid_arthritis/UncompressedCEL")
####################### NORMALISATION

# PerformING RMA normalization
normalized.data <- rma(raw.data)

# Extracting expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

gse <- getGEO("GSE13837", GSEMatrix = TRUE)

# If there are multiple series matrices, select the one you need
if (length(gse) > 1) {
  idx <- which(sapply(gse, function(x) dim(x@phenoData@data)[1]) == max(sapply(gse, function(x) dim(x@phenoData@data)[1])))
  gse_matrix <- gse[[idx]]
} else {
  gse_matrix <- gse[[1]]
}

# Now extracting feature data from the correct GEO series matrix object
feature.data <- pData(featureData(gse_matrix))

# Selecting relevant columns, assuming the 1st column is Probe ID and the 11th column is Gene Symbol

feature.data <- feature.data[, c(1, 11)]
feature.data 

# Ensuring column names in feature.data are appropriate for merging
colnames(feature.data) <- c("ID", "GeneSymbol")

# Merging the expression data with the feature data
normalized.expr <- normalized.expr %>%
  tibble::rownames_to_column(var = 'ID') %>%
  inner_join(., feature.data, by = 'ID')

# Ensuring normalized.expr is a data frame and remove any non-numeric columns
normalized.expr <- as.data.frame(normalized.expr)
numeric_columns <- sapply(normalized.expr, is.numeric)

# Keeping only numeric columns for MAD calculation
expression_data <- normalized.expr[, numeric_columns]

# Calculating the MAD for each gene
mad_values <- apply(expression_data, 1, mad, na.rm = TRUE)

# Adding the MAD values to the original data frame
normalized.expr$MAD <- mad_values

# Sorting the data by MAD in descending order
sorted_data <- normalized.expr[order(-normalized.expr$MAD), ]
# Select the top 100 genes based on MAD
top100_genes <- head(sorted_data, 100)
top100_genes


# Excluding the last two columns (MAD and GeneSymbol)
clustering_data <- top100_genes[, -seq(from = ncol(top100_genes) - 1, to = ncol(top100_genes))]
clustering_data
# exclude ID
clustering_data <- clustering_data[, -1]
clustering_data

# Specifying the file path to save the CSV file
file_path <- "~/sweta/raw_files/Rhematoid_arthritis/UncompressedCEL/clustering_data.csv"

# Saving the dataframe to a CSV file
write.csv(clustering_data, file = file_path, row.names = FALSE)
# Assuming 'clustering_data' is already loaded and it's a data frame
# First, check if 'clustering_data' is a data frame
if (!is.data.frame(clustering_data)) {
  clustering_data <- as.data.frame(clustering_data)
}

# Identifying numeric columns
numeric_columns <- sapply(clustering_data, is.numeric)

# Ensuring if numeric column is present
if (any(numeric_columns)) {
  numeric_data <- clustering_data[, numeric_columns]
} else {
  stop("No numeric columns found in the data")
}

# heatmap generation
if (exists("numeric_data")) {
  library(pheatmap)
  pheatmap(numeric_data, 
           main = "Heatmap of Top 100 Genes", 
           fontsize = 8)
}

# Converting expr_top_100_genes to a numeric matrix if needed
clustering_data <- as.matrix(clustering_data)

# Re-runing dist() to generate the distance matrix
dist_matrix <- dist(clustering_data, method = "euclidean")
dist_matrix


# Performing hierarchical clustering---- I used 5 clustering as well
hclust_result <- hclust(dist_matrix, method = "complete")
hclust_result
# Determine the number of clusters you want
num_clusters <- 10  # Adjust as needed

#DENDROGRAM generation 
cluster_assignments <- cutree(hclust_result, k = num_clusters)

# Reordering  rows and columns based on cluster assignments
reordered_data <- clustering_data[order(cluster_assignments), ]
reordered_data

# Installing  and load the pheatmap package
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)

# Assuming `clustering_data` contains the expression data for the top 100 genes
# Making sure the data is prepared correctly: rows are genes and columns are samples

# Performing  hierarchical clustering to generate  a heatmap
pheatmap(clustering_data,
         scale = "row",  # Scale the data by row
         clustering_distance_rows = "euclidean",  # Distance measure for rows
         clustering_distance_cols = "euclidean",  # Distance measure for columns
         clustering_method = "complete",  # Clustering method
         show_rownames = TRUE,  # Show or hide row names (gene names)
         show_colnames = TRUE  # Show or hide column names (sample names)
)


#####################################################


# Ensuring the row names are set correctly for clustering_data
rownames(clustering_data) <- top100_genes$ID

# Creating the annotation data frame
annotation_df <- data.frame(GeneSymbol = feature.data$GeneSymbol)
rownames(annotation_df) <- feature.data$ID
annotation_df
# Now checking and keeping only the matching row names in the annotation_df
matching_ids <- intersect(rownames(clustering_data), rownames(annotation_df))
matching_ids
annotation_df <- annotation_df[matching_ids, , drop = FALSE]
annotation_df
# If `matching_ids` is empty, there's a deeper issue with ID matching that needs to be addressed

#matching IDs
if (length(matching_ids) == 0) {
  stop("No matching gene IDs found between clustering data and annotations.")
} else {
  # Creating the heatmap with the matching annotations
  pheatmap(clustering_data[matching_ids, , drop = FALSE],  # to Ensuring only matching genes are included
           annotation_row = annotation_df,  # Adding gene annotations
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Heatmap of Top 100 Genes"
  )
}


#####################################  plotly
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}

library(plotly)
# Converting the expression matrix to a long format for plotting with plotly
clustering_data_long <- reshape2::melt(as.matrix(clustering_data))

# Creating a plotly interactive heatmap
interactive_heatmap <- plotly::plot_ly(x = clustering_data_long$Var2, # Sample names or conditions
                                       y = clustering_data_long$Var1, # Gene symbols
                                       z = clustering_data_long$value, # Expression values
                                       type = "heatmap",
                                       colorscale = "Viridis") %>% # Choose a color scale that fits your data
  plotly::layout(yaxis = list(title = "Gene Symbol", automargin = TRUE),
                 xaxis = list(title = "Samples", automargin = TRUE))

# To view the plot in an RStudio viewer or web browser
plotly::ggplotly(interactive_heatmap)

#######################  plotly with annotation

library(reshape2)
library(plotly)

# Making sure clustering_data has gene symbols as row names
rownames(clustering_data) <- top100_genes$GeneSymbol

# Melting the data for plotting with plotly
clustering_data_long <- melt(clustering_data)

# Interactive plot
interactive_heatmap <- plot_ly(x = clustering_data_long$Var2,  # These are your sample names
                               y = clustering_data_long$Var1,  # These are your gene symbols
                               z = clustering_data_long$value, # These are the expression values
                               type = "heatmap",
                               colorscale = list(c(0, "blue"), 
                                                 c(0.5, "white"), 
                                                 c(1, "red"))) %>%
  layout(yaxis = list(title = "Gene Symbol", automargin = TRUE, tickmode = "array", tickvals = clustering_data_long$Var1, ticktext = clustering_data_long$Var1),
         xaxis = list(title = "Samples", automargin = TRUE))

# Viewing the plot in RStudio or export it as HTML
ggplotly(interactive_heatmap)


# Ensuring that htmlwidgets library is present
if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets")
}
library(htmlwidgets)

# Saving html plot
saveWidget(interactive_heatmap, file = "interactive_heatmap.html")

########################################################################
library(plotly)
# First, making sure the row names are set with the gene symbols
rownames(clustering_data) <- top100_genes$GeneSymbol

# Then, melting the data for plotly
clustering_data_long <- reshape2::melt(clustering_data)

# Creating the interactive heatmap
interactive_heatmap <- plot_ly(x = clustering_data_long$Var2,  # Sample names
                               y = clustering_data_long$Var1,  # Gene symbols
                               z = clustering_data_long$value, # Expression values
                               type = "heatmap",
                               colorscale = "Viridis") %>%
  layout(yaxis = list(title = "Gene Symbol", tickfont = list(size = 8), automargin = TRUE), # Smaller font size for gene symbols
         xaxis = list(title = "Samples", automargin = TRUE))
ggplotly(interactive_heatmap)
# Saving the plot as an HTML file
htmlwidgets::saveWidget(interactive_heatmap, "interactive_heatmap.html")

# This HTML file now has a zoom feature and a smaller font size for gene symbols
#####################################################




