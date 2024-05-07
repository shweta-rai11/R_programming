library(ggplot2)
library(pheatmap)
library(reshape2)
library(arrayQualityMetrics)
library(Biobase)
library(limma)
library(MatrixGenerics)
library(RColorBrewer)
library(vsn)


#setwd("~/sweta/GSE98918")

# Load GSE42842 data set using GEOquery, here I loaded saved 
# data object from my system
load("~/sweta/GSE39428/GSE39428.RData")
show(eset)

# --------------------------------------------------------------------# 
# Examine metadata
head(pData(eset))

# Select clinico-pathological variables from sample metadata
sample_data <- pData(eset)
unique(sample_data$title)

sample_data <- data.frame( Gender = sample_data$characteristics_ch1.2,
                           state = sample_data$title,
                           treatment =sample_data$characteristics_ch1 )
head(sample_data)
tail(sample_data)


unique(sample_data$state)



# Modify the variable definitions to make them more palatable
sample_data$state[sample_data$state == "Low disease activity RA patient number 2"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 9"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 41"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 56"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 62"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 65"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 67"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 75"] <- "Low disease activity RA "
sample_data$state[sample_data$state == "Low disease activity RA patient number 89"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 99"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 104"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 113"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 138"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 147"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 154"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "Low disease activity RA patient number 155"] <- "Low disease activity RA"
sample_data$state[sample_data$state == "no Low disease activity RA patient number 19"] <- "no Low disease activity RA"
sample_data$state[sample_data$state == "no Low disease activity RA patient number 28"] <- "no Low disease activity RA"
sample_data$state[sample_data$state == "no Low disease activity RA patient number 33"] <- "no Low disease activity RA"
sample_data$state[sample_data$state == "no Low disease activity RA patient number 42"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 48"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 63"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 72"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 73"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 93"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 96"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 97"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 110"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 116"] <- "no Low disease activity RA"
sample_data$state[sample_data$state == "no Low disease activity RA patient number 117"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 118"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 127"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 137"] <- "no Low disease activity RA"

sample_data$state[sample_data$state == "no Low disease activity RA patient number 36156"] <- "no Low disease activity RA"


sample_data$state[sample_data$state == "no Low disease activity RA patient number 164"] <- "no Low disease activity RA"



sample_data$treatment[sample_data$treatment == "treatment: methotrexate /abatacept (Aba)"] <- "methotrexate /abatacept"






sample_data$Gender[sample_data$Gender == "gender: F"] <- "Female"
sample_data$Gender[sample_data$Gender == "gender: M"] <- "Male"


#sample_data$tissue[sample_data$tissue == "tissue: CD4+ memory t cell"] <- "memory_Tcell"
#sample_data$tissue[sample_data$tissue == "tissue: CD4+ naïve t cell"] <- "naive_Tcell"

#sample_data$Gender <- gsub("gender: " , "", sample_data$Gender)

#sample_data$age <- gsub("age: " , "", sample_data$age)
#class(sample_data$age) # age is still a character, we need to turn it into numbers
#sample_data$age <- as.numeric(sample_data$age)
#class(sample_data$age)

# Crosstabulate sample data
table(sample_data$state, sample_data$Gender)
# Female Male
# Control     10    5
# RA          14    4
# The samples appear balanced relative to Gender, but the majority of samples
# are female

# Convert age to categorical variable
#sort(sample_data$age)
# 45 - 80 years age range, so big difference in patient ranges
# We can separate the patients into middle-aged and old categories
#age_group<-cut(sample_data$age, 
#breaks=c(44, 58, 71), right = T)
#age_group
#sample_data$age_group <- age_group
table(sample_data$state, sample_data$Gender)
# (45,58] (58,81]
# Control      10       4
# RA           8      10
# Note that there are only 4 patients in old age control category

# -------------------------------------------------------------------------- #
# Explore array intensity distribution
# Get expression matrix
expr <- exprs(eset)
head(expr)
View(expr)


expr_df <- as.data.frame(expr)
expr_df$gene_id <- rownames(expr)  # Assuming row names are gene IDs
expr_melt <- melt(expr_df, id.vars = "id")
expr_melt

# Melt the expression matrix
expr_melt <- melt(as.data.frame(expr))
names(expr_melt) <- c("sample", "value")
head(expr_melt)


# Load the dplyr package
library(dplyr)

# Replace NA values in the "value" column with the mean
expr_melt <- expr_melt %>%
  mutate(value = ifelse(is.na(value), mean(value, na.rm = TRUE), value))

View(expr_melt)

p <- ggplot(expr_melt, aes(x=sample, y = log2(value))) +
  geom_violin() + 
  geom_boxplot(width = 0.1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
ggsave("intensities_violin-1.pdf", device = "pdf", units = "in", 
       width = 8, height = 5) 


# Check for missing values in the expression data
missing_values_exprs <- sum(is.na(exprs(eset)))

# Print the number of missing values
print(missing_values_exprs)

# Impute missing values in the expression data with the median
exprs_filled <- exprs(eset)
for (i in 1:ncol(exprs_filled)) {
  exprs_filled[is.na(exprs_filled[, i]), i] <- median(exprs_filled[, i], na.rm = TRUE)
}

# Replace the expression data in the ExpressionSet object with the imputed values
exprs(eset) <- exprs_filled

# Run array quality control analysis using arrayQualityMetrics
arrayQualityMetrics(eset,
                    outdir = "GSE68215_arrayQualityReport",
                    force = TRUE,
                    do.logtransform = TRUE,
                    intgroup = c("characteristics_ch1"),
                    spatial = TRUE,
                    reporttitle = paste("arrayQualityMetrics report for", 
                                        deparse(substitute(eset)))
)


# -------------------------------------------------------------------------- #
# Cluster the samples
# Transform gene expression values using Variance Stabilizing Transformation
# Cluster the samples
# Transform gene expression values using Variance Stabilizing Transformation
expr <- exprs(eset)
meanSdPlot(expr)
expr <- normalizeVSN(expr)
head(expr)
meanSdPlot(expr)

# Select top 500 genes with highest Mean Absolute Deviations (MAD)
mads_sorted <- sort(rowMads(expr), decreasing = TRUE)
mads_sorted <- mads_sorted[1:500]
head(mads_sorted)

# Extract gene names from mads_sorted
top_genes <- names(mads_sorted)
top_genes
str(mads_sorted)

# Extract gene names from row names of the expression data
all_genes <- rownames(exprs(eset))

# Use the top 500 MAD values to subset the gene names
top_genes <- all_genes[order(-mads_sorted)][1:500]

# Subset the expression data to include only the top 500 genes
top500 <- exprs(eset)[top_genes, ]
top500

head(sample_data)
rownames(sample_data) <- colnames(exprs(eset))

pheatmap(top500,
         scale = "none",
         clustering_method = "ward.D2",
         annotation_col = sample_data[, -3],
         show_rownames = FALSE)

# Check for missing values
any(is.na(top500))
# Check for infinite values
any(is.infinite(top500))
# Check for NA/NaN/Inf in annotation data
any(is.na(sample_data))
any(is.nan(sample_data))
any(is.infinite(sample_data))

# Check for NA values in sample_data
any(is.na(sample_data))

# Check for NaN values in each column of sample_data
apply(sample_data, 2, function(x) any(is.nan(x)))

# Check for infinite values in each column of sample_data
apply(sample_data, 2, function(x) any(is.infinite(x)))



pdf("GSE68215_heatmap_top500.pdf", width = 7, height = 6)
pheatmap(top500, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         annotation_col = sample_data[,-3], 
         show_rownames = F)
dev.off()

##################
#differential gene expression 

