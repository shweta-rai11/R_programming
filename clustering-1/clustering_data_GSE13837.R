# Download existing RA datasets

# Study 1

#Adapted Boolean Network Models for Extracellular Matrix Formation
# 	
#Wollbold J, Huber R, Pohlers D, Koczan D et al. Adapted Boolean network models for extracellular matrix formation. BMC Syst Biol 2009 Jul 21;3:77. PMID: 19622164
#Kupfer P, Guthke R, Pohlers D, Huber R et al. Batch correction of microarray data substantially improves 
#the identification of genes differentially expressed in rheumatoid arthritis and osteoarthritis. 
#BMC Med Genomics 2012 Jun 8;5:23. PMID: 22682473
#
# GEO accession number : 	GSE13837

#References used :-
# Blog post- https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
# Youtube-  https://www.youtube.com/watch?v=VMyPcUuAPyE
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5337204/
# Hierchial clustering videos by - Roger peng 
#https://www.youtube.com/watch?v=ZQYLGS7ptWM
#https://cales.arizona.edu/microarray/Jan07Workshop/lectures/Lecture%2010.pdf
#https://docs.rc.fas.harvard.edu/wp-content/uploads/2012/11/Microarray_with_R_and_bioconductor.pdf


#STEP 1: INSTALLING THE PACKAGES 

install.packages("janitor")
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  install.packages("matrixStats")
}
#STEP 2: loading the libraries
library(matrixStats)  # package to access rowSds
library(dplyr)
library(pheatmap)
library(janitor)
library(stringr)
library(dplyr)

#STEP 3: Loading the expression data
expression_data <- read.csv("GSE13837_expressions.csv", row.names = 1)
expression_data <-t(expression_data)
# Loading the metadata
metadata <- read.csv("GSE13837_phenodata.csv")

#STEP 4: Extracting the metadata- information about the sample, it is in the title column of metadata or phenodata
extract_first_word <- function(text) {
  unlist(strsplit(text, " "))[1]
}
# STEP 5:Creating a column from the 'title' ("TGF" or "TNF")
metadata <- metadata %>%
  mutate(short_title = sapply(title, extract_first_word))
#Description: The title column in metadata has longer description so i splitted it.
#STEP 6: Create a mapping from sample IDs to shortened sample names
sample_id_to_short_name <- setNames(metadata$short_title, metadata$geo_accession)
sample_id_to_short_name

#STEP 7:Changing the row names in the expression data to sample names
rownames(expression_data) <- sample_id_to_short_name[rownames(expression_data)]
View(sample_id_to_short_name)


#STEP 8: Transposing the expression data before Z-score normalization
expression_data <-t(expression_data)
normalized_expression_data <- as.data.frame(scale(expression_data, center = TRUE, scale = TRUE))

#STEP 9: Calculating MAD for each row (sample) in the normalized expression data for high variance.
mad_values_rows <- apply(normalized_expression_data, 1, function(x) {
  mad(x, na.rm = TRUE)  #missing value handling
})

#STEP 10: Sorting the samples by MAD in descending order
sorted_mad_rows <- sort(mad_values_rows, decreasing = TRUE)

#STEP 11:  Selecting only the top 500 samples with the highest MAD
top_500_samples <- names(sorted_mad_rows)[1:500]
top_500_samples

#STEP 12: normalised expression data is being subsetized.
top_500_expression_data_samples <- normalized_expression_data[top_500_samples, ]
top_500_expression_data_samples

#STEP 13:  Checking the dimensions of the subset
dim(top_500_expression_data_samples) # it returns 500

#STEP 14: Displaying the top 500 samples and their MAD values
top_500_mad_rows <- sorted_mad_rows[1:500]
top_500_mad_rows


#STEP 15: Creating dataframe FOR top_500_expression_data_samples
top_500_expression_data_samples <- as.data.frame(top_500_expression_data_samples)
top_500_expression_data_samples
dim(top_500_expression_data_samples)

#STEP 16: Saving the data frame to a CSV file for future use.
write.csv(top_500_expression_data_samples, "top_500_expression_data_samples.csv", row.names = TRUE)

############pheatmap creating

#creating pheatmap using top_500_expression_data_samples
pheatmap(top_500_expression_data_samples)
top_500_expression_data_samples
## According to the studies it is mentioned that it is recommended
#to use the scaled expression data for creating heatmaps.

# STEP 17: scaling the data
top_50_scaled <- scale(top_500_expression_data_samples)
top_50_scaled

#STEP 18: Scaling the rows
rowMeans(top_50_scaled)%>% head()
rowSds(as.matrix(top_50_scaled)) %>% head()

#checking the heatmap for entire data amd rows scaled data.
pheatmap(top_50_scaled)
pheatmap(top_500_expression_data_samples,
         scale='row')

#STEP 19: Performing hierchial clustering using Euclidean distance and linkage=complete.

my_clust <- hclust(dist(top_50_scaled, method = "euclidean"), method = "complete")  # Ensure proper linkage method
my_clust

#STEP 20: Ploting a dendrogram horizontally
as.dendrogram(my_clust)%>%
  plot(horiz=TRUE)

# STEP 21: Cutting the hierarchical clustering tree into 2 clusters only
my_gene_col <- cutree(my_clust, k = 2)
my_gene_col

#creating a dataframe 
gene_df <- data.frame(
  cluster_id = paste0("cluster", my_gene_col)
)
gene_df
#annoting the rows with the gene_df
rownames(gene_df) <- names(my_gene_col)
gene_df
print(gene_df)
#creating heatmap for rows
pheatmap(top_50_scaled,
         annotation_row = gene_df)
#for columns- removed the duplicates
colnames(top_50_scaled)
sample_col <- data.frame(
  tissue = str_remove(colnames(top_50_scaled), "\\_.")
  )
#for rows
View(sample_col)
rownames(sample_col) <- colnames(top_50_scaled)
View(sample_col)

#since it had duplicate values for the sample names so i had to remove the duplicates.
# Function to ensure unique row names
make_unique_row_names <- function(row_names) {
  data.frame(original_names = row_names) %>%
    group_by(original_names) %>%
    mutate(duplicate_count = row_number()) %>%
    ungroup() %>%
    mutate(unique_names = ifelse(duplicate_count > 1, 
                                 paste(original_names, duplicate_count, sep = "_"), 
                                 original_names)) %>%
    pull(unique_names)
}

#Dataframe with potentially duplicated row names
sample_col <- data.frame(
  tissue = gsub("_.*", "", colnames(top_50_scaled))  
)

#Ensuring the row names are unique
rownames(sample_col) <- make_unique_row_names(colnames(top_50_scaled))

# Function to ensure unique row names
make_unique_row_names <- function(row_names) {
  data.frame(original_names = row_names) %>%
    group_by(original_names) %>%
    mutate(duplicate_count = row_number()) %>%
    ungroup() %>%
    mutate(unique_names = ifelse(duplicate_count > 1, 
                                 paste(original_names, duplicate_count, sep = "_"), 
                                 original_names)) %>%
    pull(unique_names)
}
sample_col <- data.frame(
  tissue = gsub("_.*", "", colnames(top_50_scaled))  # Simplified name extraction
)

# Ensuring row names are unique
rownames(sample_col) <- make_unique_row_names(colnames(top_50_scaled))

#checking 
View(sample_col)

# Step 22: Adding sample information to the heatmap
pheatmap(top_50_scaled,
         annotation_row =gene_df,
         annotation_col = sample_col)


#cutree helps to add breaks in the rows or column for better visualization
pheatmap(top_50_scaled,
         annotation_row =gene_df,
         annotation_col = sample_col,
         cutree_rows = 2,
         cutree_cols = 2)

##################### Gene annotation
