# D. Dhaliwal, C. Sgorbati, S. Rehman 
# Install and load required packages
install.packages("naniar")
install.packages("mice")
install.packages("ggplot2")
install.packages("dplyr") # For sampling
install.packages("missForest")
install.packages("hot.deck")
install.packages("Hmisc")
install.packages("missForest")
install.packages("caret")
install.packages("Metrics")
install.packages("methods")
install.packages("tools")

# download packages for impute
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Remove the existing DESeq2 package if it exists
if ("DESeq2" %in% installed.packages()) {
  remove.packages("DESeq2")
}

# Use BiocManager to force re-install DESeq2 with dependencies
BiocManager::install("DESeq2", force = TRUE, dependencies = TRUE)

BiocManager::install("impute")
library(impute)
library(BiocManager)
library(impute)
library(caret)
library(naniar)
library(mice)
library(ggplot2)
library(dplyr)
library(missForest)
library(hot.deck)
library(Hmisc)
library(Metrics)
library(DESeq2)
library(tibble)
library(gprofiler2)

BiocManager::install('EnhancedVolcano') # package for pretty volcano plots
library(EnhancedVolcano)



## specify path 
my_path <- "/cloud/project/"

# 5% MISSING DATA
# read the data 
countdata5 <- read.table("countdata5.txt", header = T, row.names = 1,na.strings = c("NA", ""))
# print first few rows to see how the data is formatted
head(countdata5)
# transpose the data and convert back to dataframe
countdata5 <- t(countdata5)
countdata5 <- as.data.frame(countdata5)
# print total number of missing values
total_missing <- sum(is.na(countdata5))
print(total_missing)
# remove columns with only zeroes
countdata5_clean <- countdata5[ , colSums(countdata5 != 0,na.rm=T)>0] 
# check how many NA values there are
clean_missing <- sum(is.na(countdata5_clean))
print(clean_missing)
# calculate the proportion of the data that is missing values to check it is still 5%
prop_miss(countdata5_clean)
# calculate the sum of each column
colSums(countdata5_clean)
# create upset plot
gg_miss_upset(countdata5_clean)
# print the summary of missing values
miss_summary_5 <- miss_var_summary(countdata5_clean)
print(miss_summary_5)
# other visualisations
gg_miss_var(countdata5_clean)
gg_miss_case(countdata5_clean)
gg_miss_span(countdata5_clean,
             span_every = 1)

#kNN imputation
# Convert the count data to a matrix.
# This is necessary because the impute.knn function requires a matrix input.
matCountdata5 <- as.matrix(countdata5_clean)
# Perform k-nearest neighbors imputation on the count data matrix.
# The impute.knn function imputes missing values using the k-nearest neighbors algorithm.
# Here, k = 10 means that the 10 nearest neighbors are used for imputation.
knn_test_5 <- impute.knn(matCountdata5, k = 10)
# Convert the resulting imputed matrix back to a data frame.
knn_dataframe_5 <- as.data.frame(knn_test_5$data)
# Display the first few rows of the imputed data frame.
head(knn_dataframe_5)
# Check the number of remaining NA values in the imputed data frame.
# This should be zero if the imputation was successful.
sum(is.na(knn_dataframe_5))
# Save the imputed data frame as a tab-separated text file.
# This allows for easy sharing and future use of the imputed data.
write.table(knn_dataframe_5, paste0(my_path, "knn_countdata5.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
sum(is.na(knn_dataframe_5))


#RandomForest
# Perform random forest imputation on the count data matrix.
# The missForest function uses the random forest algorithm to estimate and fill in missing values.
rf_5_test <- missForest(matCountdata5)
# Convert the resulting imputed matrix back to a data frame.
rf_5_dataframe <- as.data.frame(rf_5_test$ximp)
# Display the first few rows of the imputed data frame.
head(rf_5_dataframe)
# Check the number of remaining NA values in the imputed data frame.
# This should be zero if the imputation was successful.
sum(is.na(rf_5_dataframe))
# Save the imputed data frame as a tab-separated text file.
# This allows for easy sharing and future use of the imputed data.
write.table(rf_5_dataframe, paste0(my_path, "rf5_countdata.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
# save as a spreadsheet
write.table(rf_5_dataframe, paste0(my_path, "rf_countdata5_randomforest.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
sum(is.na(rf_5_dataframe))


#mean imputation
mean_imputation_5 <- mice(matCountdata5, method = 'mean', m = 5, maxit = 10)
mean_imputation_5 <- complete(mean_imputation_5, 1)
mean_imputation_5 <- as.data.frame(mean_imputation_5)
summary(mean_imputation_5)
# save as a spreadsheet
write.table(mean_imputation_5, paste0(my_path, "mean_countdata5.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)


# DESeq2

# for kNN
knn5_countdata <- read.table(paste0(my_path, "knn_countdata5.txt"), header = TRUE, row.names = 1)
coldata <- read.csv(paste0(my_path, "metadata.csv"))
# Define groups based on latitude
coldata$group <- "North"
coldata$group[coldata$Latitude < 35] <- "South"
# Keep only the distant populations
coldata <- rbind(coldata[coldata$Latitude < 35, ], coldata[coldata$Latitude > 40, ])
# Subset knn5_countdata to match coldata
knn5_countdata <- knn5_countdata[, colnames(knn5_countdata) %in% coldata$RNAseq_Run_NCBI]
# Ensure the columns of knn5_countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(knn5_countdata), coldata$RNAseq_Run_NCBI), ]
# Check if the columns match after subsetting
if (!all(colnames(knn5_countdata) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}
# Remove rows with NA values
knn5_countdata <- na.omit(knn5_countdata)
# Round the count data to ensure all values are integers
knn5_countdata <- round(knn5_countdata)
# Create a DESeq object
knn_dds_5 <- DESeqDataSetFromMatrix(countData = knn5_countdata,
                                   colData = coldata,
                                   design = ~ Sex + group)
# Proceed with the differential expression analysis
knn_dds_5 <- DESeq(knn_dds_5)
# Results
res5knn <- results(knn_dds_5)
results_5knn <- as.data.frame(res5knn)
results_5knn$Ensembl_id <- row.names(results_5knn)
results_5knn <- results_5knn[order(results_5knn$padj), ]
# Sanity check: Number of genes before the merge
num_genes_beforemerge_5knn <- nrow(results_5knn)
cat("Number of genes before merge:", num_genes_beforemerge_5knn, "\n")
# Convert Ensembl names to gene names (if they are known)
results_genes_5knn <- gconvert(row.names(res5knn), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)
# Ensure unique entries in results_genes_5knn
results_genes_5knn <- results_genes_5knn[!duplicated(results_genes_5knn$input), ]
# Add the gene names
results_5knn <- merge(results_5knn,
                      results_genes_5knn[, c("input", "target", "name", "description")],
                      by.x = "Ensembl_id", by.y = "input", all.x = TRUE)
# Sanity check: Number of genes after the merge
num_genes_beforemerge_5knn <- nrow(results_5knn)
cat("Number of genes after merge:", num_genes_beforemerge_5knn, "\n")
# Ensure the 'Name' column is properly populated
results_5knn$Name <- ifelse(is.na(results_5knn$name), results_5knn$Ensembl_id, results_5knn$name)
# Check for missing values in the 'Name' column
num_na_names_5knn <- sum(is.na(results_5knn$Name))
cat("Number of missing names after merge:", num_na_names_5knn, "\n")
# Subset the results to keep only significant genes
significant_5knngenes <- results_5knn[results_5knn$padj < 0.05 & !is.na(results_5knn$padj), ]
# Ensure significant_5knngenes is not empty before plotting
if (nrow(significant_5knngenes) > 0) {
  # Volcano plot
  EnhancedVolcano(significant_5knngenes,
                  lab = significant_5knngenes$Name,
                  x = 'log2FoldChange',
                  y = 'padj')
} else {
  message("No significant genes found.")
}
# Check the significant genes
head(significant_5knngenes)


# for MissForest
# Read in the count data (rf 5%)
rf5_countdata <- read.table(paste0(my_path, "rf_countdata5_randomforest.txt"), header = TRUE, row.names = 1)
# Subset rf5_countdata to match coldata
rf5_countdata <- rf5_countdata[, colnames(rf5_countdata) %in% coldata$RNAseq_Run_NCBI]
# Ensure the columns of rf5_countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(rf5_countdata), coldata$RNAseq_Run_NCBI), ]
# Check if the columns match after subsetting
if (!all(colnames(rf5_countdata) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}
# Remove rows with NA values
rf5_countdata <- na.omit(rf5_countdata)
# Round the count data to ensure all values are integers
rf5_countdata <- round(rf5_countdata)
# Create a DESeq object
rf5_dds <- DESeqDataSetFromMatrix(countData = rf5_countdata,
                                  colData = coldata,
                                  design = ~ Sex + group)
# Proceed with the differential expression analysis
rf5_dds <- DESeq(rf5_dds)
# Results
res_rf5 <- results(rf5_dds)
results_rf5 <- as.data.frame(res_rf5)
results_rf5$Ensembl_id <- row.names(results_rf5)
results_rf5 <- results_rf5[order(results_rf5$padj), ]
# Sanity check: Number of genes before the merge
num_genes_beforemerge_rf5 <- nrow(results_rf5)
cat("Number of genes before merge:", num_genes_beforemerge_rf5, "\n")
# Convert Ensembl names to gene names (if they are known)
results_genes_rf5 <- gconvert(row.names(res_rf5), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)
# Ensure unique entries in results_genes_rf5
results_genes_rf5 <- results_genes_rf5[!duplicated(results_genes_rf5$input), ]
# Add the gene names
results_rf5 <- merge(results_rf5,
                     results_genes_rf5[, c("input", "target", "name", "description")],
                     by.x = "Ensembl_id", by.y = "input", all.x = TRUE)
# Sanity check: Number of genes after the merge
num_genes_aftermerge_rf5 <- nrow(results_rf5)
cat("Number of genes after merge:", num_genes_aftermerge_rf5, "\n")
# Ensure the 'Name' column is properly populated
results_rf5$Name <- ifelse(is.na(results_rf5$name), results_rf5$Ensembl_id, results_rf5$name)
# Check for missing values in the 'Name' column
num_na_names_rf5 <- sum(is.na(results_rf5$Name))
cat("Number of missing names after merge:", num_na_names_rf5, "\n")
# Subset the results to keep only significant genes
significant_rf5_genes <- results_rf5[results_rf5$padj < 0.05 & !is.na(results_rf5$padj), ]
# Ensure significant_rf5_genes is not empty before plotting
if (nrow(significant_rf5_genes) > 0) {
  # Remove rows with NA values in the columns used for plotting
  significant_rf5_genes <- significant_rf5_genes[!is.na(significant_rf5_genes$log2FoldChange) & !is.na(significant_rf5_genes$padj), ]
  # Check the range of log2FoldChange and padj values
  cat("Range of log2FoldChange:", range(significant_rf5_genes$log2FoldChange, na.rm = TRUE), "\n")
  cat("Range of padj:", range(significant_rf5_genes$padj, na.rm = TRUE), "\n")
  # Proceed with the volcano plot if the cleaned dataset is not empty
  if (nrow(significant_rf5_genes) > 0) {
    # Volcano plot
    EnhancedVolcano(significant_rf5_genes,
                    lab = significant_rf5_genes$Name,
                    x = 'log2FoldChange',
                    y = 'padj')
  } else {
    message("No significant genes found after cleaning.")
  }
} else {
  message("No significant genes found.")
}
# Check the significant genes
head(significant_rf5_genes)



# for mean
# Read in the count data (mean 5%) 
mean5_countdata <- read.table(paste0(my_path, "mean_countdata5.txt"), header = TRUE, row.names = 1)
# Subset mean5_countdata to match coldata
mean5_countdata <- mean5_countdata[, colnames(mean5_countdata) %in% coldata$RNAseq_Run_NCBI]
# Ensure the columns of mean5_countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(mean5_countdata), coldata$RNAseq_Run_NCBI), ]
# Check if the columns match after subsetting
if (!all(colnames(mean5_countdata) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}
# Remove rows with NA values
mean5_countdata <- na.omit(mean5_countdata)
# Round the count data to ensure all values are integers
mean5_countdata <- round(mean5_countdata)
# Create a DESeq object
mean5_dds <- DESeqDataSetFromMatrix(countData = mean5_countdata,
                                    colData = coldata,
                                    design = ~ Sex + group)
# Proceed with the differential expression analysis
mean5_dds <- DESeq(mean5_dds)
# Results
res_mean5 <- results(mean5_dds)
results_mean5 <- as.data.frame(res_mean5)
results_mean5$Ensembl_id <- row.names(results_mean5)
results_mean5 <- results_mean5[order(results_mean5$padj), ]
# Sanity check: Number of genes before the merge
num_genes_beforemerge_mean5 <- nrow(results_mean5)
cat("Number of genes before merge:", num_genes_beforemerge_mean5, "\n")
# Convert Ensembl names to gene names (if they are known)
results_genes_mean5 <- gconvert(row.names(res_mean5), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)
# Check for missing mappings in gconvert output
missing_genes_mean5 <- setdiff(row.names(res_mean5), results_genes_mean5$input)
cat("Number of missing genes:", length(missing_genes_mean5), "\n")
# Ensure unique entries in results_genes_mean5
results_genes_mean5 <- results_genes_mean5[!duplicated(results_genes_mean5$input), ]
# Add the gene names
results_mean5 <- merge(results_mean5,
                       results_genes_mean5[, c("input", "target", "name", "description")],
                       by.x = "Ensembl_id", by.y = "input", all.x = TRUE)
# Sanity check: Number of genes after the merge
num_genes_aftermerge_mean5 <- nrow(results_mean5)
cat("Number of genes after merge:", num_genes_aftermerge_mean5, "\n")
# Ensure the 'Name' column is properly populated
results_mean5$Name <- ifelse(is.na(results_mean5$name), results_mean5$Ensembl_id, results_mean5$name)
# Check for missing values in the 'Name' column
num_na_names_mean5 <- sum(is.na(results_mean5$Name))
cat("Number of missing names after merge:", num_na_names_mean5, "\n")
# Subset the results to keep only significant genes
significant_mean5_genes <- results_mean5[results_mean5$padj < 0.05 & !is.na(results_mean5$padj), ]
# Ensure significant_mean5_genes is not empty before plotting
if (nrow(significant_mean5_genes) > 0) {
  # Remove rows with NA values in the columns used for plotting
  significant_mean5_genes <- significant_mean5_genes[!is.na(significant_mean5_genes$log2FoldChange) & !is.na(significant_mean5_genes$padj), ]
  # Check the range of log2FoldChange and padj values
  cat("Range of log2FoldChange:", range(significant_mean5_genes$log2FoldChange, na.rm = TRUE), "\n")
  cat("Range of padj:", range(significant_mean5_genes$padj, na.rm = TRUE), "\n")
  # Proceed with the volcano plot if the cleaned dataset is not empty
  if (nrow(significant_mean5_genes) > 0) {
    # Volcano plot
    EnhancedVolcano(significant_mean5_genes,
                    lab = significant_mean5_genes$Name,
                    x = 'log2FoldChange',
                    y = 'padj')
  } else {
    message("No significant genes found after cleaning.")
  }
} else {
  message("No significant genes found.")
}

# Check the significant genes
head(significant_mean5_genes)

# EVALUATION
# evaluate knn
# create function to evaluate the knn
knn_imputation_evaluation <- function(file_path, k = 10, na_fraction = 0.3, seed1 = 123, seed2 = 456) {
  # Helper function for KNN imputation
  perform_knn_imputation <- function(data, k) {
    imputed <- impute.knn(as.matrix(data), k = k)$data
    return(as.data.frame(imputed))
  }
  
  # Load the dataset
  countdata_imputedKNN <- read.table(file_path, header = TRUE, row.names = 1)
  
  # Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed1) # For reproducibility
  total_samples <- nrow(countdata_imputedKNN)
  train_size <- floor(0.7 * total_samples)
  train_indices <- sample(seq_len(total_samples), size = train_size)
  train_data <- countdata_imputedKNN[train_indices, ]
  test_data <- countdata_imputedKNN[-train_indices, ]
  
  # Introduce artificial NAs randomly in the testing set
  set.seed(seed2) # For reproducibility
  test_data_with_na <- test_data
  num_rows <- nrow(test_data_with_na)
  num_cols <- ncol(test_data_with_na)
  num_values <- prod(dim(test_data_with_na))
  num_na <- round(num_values * na_fraction)
  
  # Calculate the maximum number of NAs allowed per column and row (80% of each)
  max_na_per_col <- floor(0.8 * num_rows)
  max_na_per_row <- floor(0.8 * num_cols)
  
  # Function to distribute NAs ensuring no row/column exceeds 80% NAs
  distribute_nas <- function(test_data_with_na, num_na, max_na_per_col, max_na_per_row) {
    na_positions <- array(FALSE, dim(test_data_with_na))
    current_na_count <- 0
    while (current_na_count < num_na) {
      # Randomly choose a position
      row_idx <- sample(num_rows, 1)
      col_idx <- sample(num_cols, 1)
      
      # Check if adding an NA here would exceed the 80% limit
      if (sum(na_positions[row_idx, ]) < max_na_per_row &&
          sum(na_positions[, col_idx]) < max_na_per_col &&
          !na_positions[row_idx, col_idx]) {
        na_positions[row_idx, col_idx] <- TRUE
        current_na_count <- current_na_count + 1
      }
    }
    return(na_positions)
  }
  
  # Distribute NAs according to the constraint
  na_positions <- distribute_nas(test_data_with_na, num_na, max_na_per_col, max_na_per_row)
  test_data_with_na[na_positions] <- NA
  
  # Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.")
  }
  
  # Re-impute the missing values in the testing set using kNN
  test_data_reimputed <- perform_knn_imputation(test_data_with_na, k = k)
  
  # Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.")
  }
  
  # Extract true and imputed values for comparison
  true_values <- test_data[is.na(test_data_with_na)]
  imputed_values <- test_data_reimputed[is.na(test_data_with_na)]
  
  # Ensure there are no NAs in the extracted true and imputed values
  if (any(is.na(true_values)) || any(is.na(imputed_values))) {
    stop("NAs present in the true or imputed values for comparison.")
  }
  
  # Calculate MAE and RMSE between the imputed values and the true values in the test set
  mae_value <- mae(true_values, imputed_values)
  rmse_value <- rmse(true_values, imputed_values)
  cat("Mean Absolute Error (MAE):", mae_value, "\n")
  cat("Root Mean Squared Error (RMSE):", rmse_value, "\n")
  
  # Save the datasets
  write.table(train_data, file = "train_knn.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_with_na, file = "test_knn.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_reimputed, file = "knn_countdata_test_reimputed.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Plot the distribution of true values and imputed values
  true_vs_imputed <- data.frame(True = true_values, Imputed = imputed_values)
  
  # Plot a scatter plot and a line y=x for reference
  p1 <- ggplot(true_vs_imputed, aes(x = True, y = Imputed)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    labs(title = "True vs Imputed Values", x = "True Values", y = "Imputed Values") +
    theme_minimal()
  
  # Plot the error distribution
  errors <- data.frame(Error = true_values - imputed_values)
  
  p2 <- ggplot(errors, aes(x = Error)) +
    geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Error Distribution (True - Imputed)", x = "Error", y = "Frequency") +
    theme_minimal()
  
  print(p1)
  print(p2)
}

      
knn_imputation_evaluation(paste0(my_path, "knn_countdata5.txt"))

# evaluate missforest
# create function to evaluate the missforest
missforest_evaluation_function <- function(data, na_fraction = 0.3, seed_split = 123, seed_na = 456) {
  # Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed_split) # For reproducibility
  total_samples <- nrow(data)
  train_size <- floor(0.7 * total_samples)
  
  train_indices <- sample(seq_len(total_samples), size = train_size)
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  # Introduce artificial NAs randomly in the testing set
  set.seed(seed_na) # For reproducibility
  test_data_with_na <- test_data
  num_values <- prod(dim(test_data_with_na))
  num_na <- round(num_values * na_fraction)
  na_indices <- arrayInd(sample(num_values, num_na), dim(test_data_with_na))
  test_data_with_na[na_indices] <- NA
  
  # Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.")
  }
  
  # Re-impute the missing values in the testing set using missForest
  imputed_data <- missForest(test_data_with_na)
  
  # Extract the reimputed data
  test_data_reimputed <- imputed_data$ximp
  
  # Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.")
  }
  
  # Extract true and imputed values for comparison
  true_values <- test_data[na_indices]
  imputed_values <- test_data_reimputed[na_indices]
  
  # Ensure there are no NAs in the extracted true and imputed values
  if (any(is.na(true_values)) || any(is.na(imputed_values))) {
    stop("NAs present in the true or imputed values for comparison.")
  }
  
  # Calculate MAE and RMSE between the imputed values and the true values in the test set
  mae_value <- mae(true_values, imputed_values)
  rmse_value <- rmse(true_values, imputed_values)
  
  cat("Mean Absolute Error (MAE):", mae_value, "\n")
  cat("Root Mean Squared Error (RMSE):", rmse_value, "\n")
  
  # Save the datasets
  write.table(train_data, file = paste0(my_path, "train_mf_30.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_with_na, file = paste0(my_path, "mf_countdata_test_with_na.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_reimputed, file = paste0(my_path, "mf_countdata_test_reimputed.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Plot the distribution of true values and imputed values
  true_vs_imputed <- data.frame(True = true_values, Imputed = imputed_values)
  
  # Plot a scatter plot and a line y=x for reference
  p1 <- ggplot(true_vs_imputed, aes(x = True, y = Imputed)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    labs(title = "True vs Imputed Values", x = "True Values", y = "Imputed Values") +
    theme_minimal()
  
  # Plot the error distribution
  errors <- data.frame(Error = true_values - imputed_values)
  
  p2 <- ggplot(errors, aes(x = Error)) +
    geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Error Distribution (True - Imputed)", x = "Error", y = "Frequency") +
    theme_minimal()
  
  print(p1)
  print(p2)
}

missforest_evaluation_function(rf_5_dataframe)

# mean evaluation
mean_imputation_evaluation <- function(file_path, na_fraction = 0.3, seed1 = 123, seed2 = 456) {
  perform_mean_imputation <- function(data) {
    imputed <- as.data.frame(lapply(data, function(x) impute(x, fun = mean)))
    return(imputed)
  }
  
  # Load the dataset
  countdata10_imputedMean <- read.table(paste0(my_path, "mean_imp.txt"), header = TRUE, row.names = 1)
  
  # Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed1) # For reproducibility
  total_samples <- nrow(countdata10_imputedMean)
  train_size <- floor(0.7 * total_samples)
  
  train_indices <- sample(seq_len(total_samples), size = train_size)
  train_data <- countdata10_imputedMean[train_indices, ]
  test_data <-  countdata10_imputedMean[-train_indices, ]
  
  # Introduce artificial NAs randomly in the testing set
  set.seed(seed2) # For reproducibility
  test_data_with_na <- test_data
  num_values <- prod(dim(test_data_with_na))
  num_na <- round(num_values * na_fraction)
  na_indices <- arrayInd(sample(num_values, num_na), dim(test_data_with_na))
  test_data_with_na[na_indices] <- NA
  
  # Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.")
  }
  
  # Impute missing values using mean imputation
  perform_mean_imputation <- function(data) {
    imputed <- as.data.frame(lapply(data, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))
    return(imputed)
  }
  
  test_data_reimputed <- perform_mean_imputation(test_data_with_na)
  
  # Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.")
  }
  
  # Extract true and imputed values for comparison
  true_values <- test_data[na_indices]
  imputed_values <- test_data_reimputed[na_indices]
  
  # Ensure there are no NAs in the extracted true and imputed values
  if (any(is.na(true_values)) || any(is.na(imputed_values))) {
    stop("NAs present in the true or imputed values for comparison.")
  }
  
  # Calculate MAE and RMSE between the imputed values and the true values in the test set
  mae_value <- mean(abs(true_values - imputed_values))
  rmse_value <- sqrt(mean((true_values - imputed_values)^2))
  
  cat("Mean Absolute Error (MAE):", mae_value, "\n")
  cat("Root Mean Squared Error (RMSE):", rmse_value, "\n")
  
  # Save the datasets
  write.table(train_data, file = paste0(my_path, "train_mean10_imp.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_with_na, file = paste0(my_path, "test_mean10_imp.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_reimputed, file = paste0(my_path, "mean10_countdata_test_reimputed.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Plot the distribution of true values and imputed values
  true_vs_imputed <- data.frame(True = true_values, Imputed = imputed_values)
  
  # Plot a scatter plot and a line y=x for reference
  p1 <- ggplot(true_vs_imputed, aes(x = True, y = Imputed)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    labs(title = "True vs Imputed Values", x = "True Values", y = "Imputed Values") +
    theme_minimal()
  
  # Plot the error distribution
  errors <- data.frame(Error = true_values - imputed_values)
  
  p2 <- ggplot(errors, aes(x = Error)) +
    geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Error Distribution (True - Imputed)", x = "Error", y = "Frequency") +
    theme_minimal()
  
  print(p1)
  print(p2)
}

# use on the mean imputed dataset 
mean_imputation_evaluation_mean5(paste0(my_path, "mean_countdata5.txt"))

# 10% MISSING DATA
# Load the dataset
countdata10 <- read.table("countdata10.txt", header = T, row.names = 1)
# transpose dataset
countdata <- as.data.frame(t(countdata10))
# remove the 0 
countdata_clean_t = countdata[ , colSums(countdata != 0, na.rm = TRUE) > 0]
# transpose back for kNN and missForest
countdata_10_clean = as.data.frame(t(countdata_clean_t))
# Visualize missing data patterns
gg_miss_upset(countdata_clean_t)
# Generate a summary of the missing values
miss_summary_10 <- miss_var_summary(countdata_clean_t)
print(miss_summary_10)
# transform dataframe into matrix (use non-transposed data)
matrixcount <- as.matrix(countdata_10_clean)
# impute with kNN
test_knn_10 <- impute.knn(data=matrixcount, k = 10)
# check for missing values
sum(is.na(test_knn_10))
# Extract the imputed dataframe from the list
knn10_imp <-  as.data.frame(test_knn_10$data)  
# save as text file
write.table(knn10_imp, "knn10_imp.txt", row.names =TRUE, col.names = TRUE, quote = FALSE)


# missforest
# impute with random forest
test_missforest_10 <- missForest(matrixcount)
# check for missing values
sum(is.na(test_missforest_10))
# Extract the imputed dataframe from the list
mf10_imp <-  as.data.frame(test_missforest_10$ximp)
# save as text file
write.table(mf10_imp, "mf10_imp.txt", row.names =TRUE, col.names = TRUE, quote = FALSE)

# impute with mean 
# Impute missing values with mean
mean_imp_10 <- as.data.frame(lapply(countdata_clean_t, function(x) impute(x, fun = mean)))
# check for missing values
sum(is.na(mean_imp_10))
# save as text file
write.table(mean_imp_10, "mean_imp_10.txt", row.names =TRUE, col.names = TRUE, quote = FALSE)

# DESeq2 
## DESeq2 knn
# Read in the count data and metadata
countdata_knn10 <- read.table(paste0(my_path, "knn10_imp.txt"), header = TRUE, row.names = 1)

# Subset countdata to match coldata
countdata_knn10 <- countdata_knn10[, colnames(countdata_knn10) %in% coldata$RNAseq_Run_NCBI]

# Ensure the columns of countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(countdata_knn10), coldata$RNAseq_Run_NCBI), ]

# Check if the columns match after subsetting
if (!all(colnames(countdata_knn10) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}

# Remove rows with NA values
countdata_knn10 <- na.omit(countdata_knn10)

# Round the count data to ensure all values are integers
countdata_knn10 <- round(countdata_knn10)

# Create a DESeq object
dds_knn10 <- DESeqDataSetFromMatrix(countData = countdata_knn10,
                                    colData = coldata,
                                    design = ~ Sex + group)
head(dds_knn10)

# Proceed with the differential expression analysis
dds_knn10 <- DESeq(dds_knn10)

# Results
res_knn10 <- results(dds_knn10)
results_df_knn10 <- as.data.frame(res_knn10)
results_df_knn10$Ensembl_id <- row.names(results_df_knn10)
results_df_knn10 <- results_df_knn10[order(results_df_knn10$padj), ]

sum(is.na(results_df_knn10))

# Sanity check: Number of genes before the merge
num_genes_before_merge <- nrow(results_df_knn10)
cat("Number of genes before merge:", num_genes_before_merge, "\n")

# Convert Ensembl names to gene names (if they are known)
results_genes_knn10 <- gconvert(row.names(res_knn10), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)

# Ensure unique entries in results_genes
results_genes_knn10 <- results_genes_knn10[!duplicated(results_genes_knn10$input), ]

# Add the gene names
results_df_knn10 <- merge(results_df_knn10,
                          results_genes_knn10[, c("input", "target", "name", "description")],
                          by.x = "Ensembl_id", by.y = "input", all.x = TRUE)

# Sanity check: Number of genes after the merge
num_genes_after_merge <- nrow(results_df_knn10)
cat("Number of genes after merge:", num_genes_after_merge, "\n")

# Ensure the 'Name' column is properly populated
results_df_knn10$Name <- ifelse(is.na(results_df_knn10$name), results_df_knn10$Ensembl_id, results_df_knn10$name)

# Check for missing values in the 'Name' column
num_na_names <- sum(is.na(results_df_knn10$Name))
cat("Number of missing names after merge:", num_na_names, "\n")

# Subset the results to keep only significant genes
significant_genes_knn10 <- results_df_knn10[results_df_knn10$padj < 0.05 & !is.na(results_df_knn10$padj), ]
head(significant_genes_knn10)

# Ensure significant_genes is not empty before plotting
if (nrow(significant_genes_knn10) > 0) {
  # Volcano plot
  EnhancedVolcano(significant_genes_knn10,
                  lab = significant_genes_knn10$Name,
                  x = 'log2FoldChange',
                  y = 'padj')
} else {
  message("No significant genes found.")
}

# save results
write.csv(results_df, file = "results_DE-knn10%.csv", row.names = F)

## DESeq2 mf 
# Read in the count data and metadata
countdata_mf10 <- read.table(paste0(my_path, "mf10_imp.txt"), header = TRUE, row.names = 1)

# Subset countdata to match coldata
countdata_mf10 <- countdata_mf10[, colnames(countdata_mf10) %in% coldata$RNAseq_Run_NCBI]

# Ensure the columns of countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(countdata_mf10), coldata$RNAseq_Run_NCBI), ]

# Check if the columns match after subsetting
if (!all(colnames(countdata_mf10) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}

# Remove rows with NA values
countdata_mf10 <- na.omit(countdata_mf10)

# Round the count data to ensure all values are integers
countdata_mf10 <- round(countdata_mf10)

# Create a DESeq object
dds_mf10 <- DESeqDataSetFromMatrix(countData = countdata_mf10,
                                   colData = coldata,
                                   design = ~ Sex + group)

# Proceed with the differential expression analysis
dds_mf10 <- DESeq(dds_mf10)

# Results
res_mf10 <- results(dds_mf10)
results_df_mf10 <- as.data.frame(res_mf10)
results_df_mf10$Ensembl_id <- row.names(results_df_mf10)
results_df_mf10 <- results_df_mf10[order(results_df_mf10$padj), ]

# Sanity check: Number of genes before the merge
num_genes_before_merge <- nrow(results_df_mf10)
cat("Number of genes before merge:", num_genes_before_merge, "\n")

# Convert Ensembl names to gene names (if they are known)
results_genes_mf10 <- gconvert(row.names(res_mf10), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)

# Ensure unique entries in results_genes
results_genes_mf10 <- results_genes_mf10[!duplicated(results_genes_mf10$input), ]

# Add the gene names
results_df_mf10 <- merge(results_df_mf10,
                         results_genes_mf10[, c("input", "target", "name", "description")],
                         by.x = "Ensembl_id", by.y = "input", all.x = TRUE)

# Sanity check: Number of genes after the merge
num_genes_after_merge <- nrow(results_df_mf10)
cat("Number of genes after merge:", num_genes_after_merge, "\n")

# Ensure the 'Name' column is properly populated
results_df_mf10$Name <- ifelse(is.na(results_df_mf10$name), results_df_mf10$Ensembl_id, results_df_mf10$name)

# Check for missing values in the 'Name' column
num_na_names <- sum(is.na(results_df_mf10$Name))
cat("Number of missing names after merge:", num_na_names, "\n")

# Subset the results to keep only significant genes
significant_genes_mf10 <- results_df_mf10[results_df_mf10$padj < 0.05 & !is.na(results_df_mf10$padj), ]
head(significant_genes_mf10)

# Ensure significant_genes is not empty before plotting
if (nrow(significant_genes_mf10) > 0) {
  # Volcano plot
  EnhancedVolcano(significant_genes_mf10,
                  lab = significant_genes_mf10$Name,
                  x = 'log2FoldChange',
                  y = 'padj')
} else {
  message("No significant genes found.")
}

# save results
write.csv(results_df, file = "results_DE-mf10%.csv", row.names = F)

## DESeq2 for mean
# Read in the count data and metadata
countdata_mean10 <- read.table(paste0(my_path, "mean_imp.txt"), header = TRUE, row.names = 1)

# Subset countdata to match coldata
countdata_mean10 <- countdata_mean10[, colnames(countdata_mean10) %in% coldata$RNAseq_Run_NCBI]

# Ensure the columns of countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(countdata_mean10), coldata$RNAseq_Run_NCBI), ]

# Check if the columns match after subsetting
if (!all(colnames(countdata_mean10) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}

# Remove rows with NA values
countdata_mean10 <- na.omit(countdata_mean10)

# Round the count data to ensure all values are integers
countdata_mean10 <- round(countdata_mean10)

# Create a DESeq object
dds_mean10 <- DESeqDataSetFromMatrix(countData = countdata_mean10,
                                     colData = coldata,
                                     design = ~ Sex + group)

# Proceed with the differential expression analysis
dds_mean10  <- DESeq(dds_mean10)

# Results
res_mean10  <- results(dds_mean10)
results_df_mean10  <- as.data.frame(res_mean10)
results_df_mean10$Ensembl_id <- row.names(results_df_mean10)
results_df_mean10 <- results_df_mean10[order(results_df_mean10$padj), ]

# Sanity check: Number of genes before the merge
num_genes_before_merge <- nrow(results_df_mean10)
cat("Number of genes before merge:", num_genes_before_merge, "\n")

# Convert Ensembl names to gene names (if they are known)
results_genes_mean10 <- gconvert(row.names(res_mean10), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)

# Ensure unique entries in results_genes
results_genes_mean10 <- results_genes_mean10[!duplicated(results_genes_mean10$input), ]

# Add the gene names
results_df_mean10 <- merge(results_df_mean10,
                           results_genes_mean10[, c("input", "target", "name", "description")],
                           by.x = "Ensembl_id", by.y = "input", all.x = TRUE)

# Sanity check: Number of genes after the merge
num_genes_after_merge <- nrow(results_df_mean10)
cat("Number of genes after merge:", num_genes_after_merge, "\n")

# Ensure the 'Name' column is properly populated
results_df_mean10$Name <- ifelse(is.na(results_df_mean10$name), results_df_mean10$Ensembl_id, results_df_mean10$name)

# Check for missing values in the 'Name' column
num_na_names <- sum(is.na(results_df_mean10$Name))
cat("Number of missing names after merge:", num_na_names, "\n")

# Subset the results to keep only significant genes
significant_genes_mean10 <- results_df_mean10[results_df_mean10$padj < 0.05 & !is.na(results_df_mean10$padj), ]

# Ensure significant_genes is not empty before plotting
if (nrow(significant_genes_mean10) > 0) {
  # Volcano plot
  EnhancedVolcano(significant_genes_mean10,
                  lab = significant_genes_mean10$Name,
                  x = 'log2FoldChange',
                  y = 'padj')
} else {
  message("No significant genes found.")
}

# save results
write.csv(results_df, file = "results_DE-mean10%.csv", row.names = F)

# evaluation 
# kNN
knn_imputation_evaluation(paste0(my_path, "knn10_imp.txt"))

# missforest 
missforest_evaluation_function(mf10_imp)

# mean 
mean_imputation_evaluation(paste0(my_path, "mean_imp.txt"))

# 30% MISSING DATA

# determine what the type of missing data is and visualise
# Load the data
countdata30 <- read.table(paste0(my_path, "countdata30.txt"), header = TRUE, row.names = 1)
# Convert all columns to numeric
countdata30 <- as.data.frame(lapply(countdata30, function(x) as.numeric(as.character(x))))
# convert into dataframe
countdata30 <- as.data.frame(countdata30)
# check it is dataframe and print first 6 rows
is.data.frame(countdata30)
head(countdata30)
# Visualize the pattern of missing data
gg_miss_var(countdata30) + labs(title = "Missing Data Pattern by Mice in Countdata 30%")
# Additional visualization options
gg_miss_upset(countdata30)
gg_miss_case(countdata30)
# Summarize missing data by variables (which are now mice)
miss_summary_30 <- miss_var_summary(countdata30)
print(miss_summary_30)


# knn imputation
# check number of NA values 
sum(is.na(countdata30))
# create function to remove columns with only zero values
remove_zero_rows <- function(data) {
  # Identify columns where all values are zero
  non_zero_rows <- rowSums(data != 0, na.rm = TRUE) > 0
  # Subset the data to keep only columns with non-zero values
  cleaned_data <- data[non_zero_rows, ]
  # Print the number of columns removed
  num_removed <- sum(!non_zero_rows)
  print(paste("Number of rows removed due to only zero values:", num_removed))
  return(cleaned_data)
}
# apply the function to countdata30 to remove only zero columns from function
countdata30_clean <- remove_zero_rows(countdata30)
# check the number of NAs after removing columns of zeroes
sum(is.na(countdata30_clean))
# calculate the proportion of the data that is missing values to check it is still 30%
prop_miss(countdata30_clean)
# convert clean dataset to a matrix 
mat_countdata30_clean <- as.matrix(countdata30_clean)
# impute for the clean dataset as matrix 
clean_kNN_imputed_30 <- impute.knn(mat_countdata30_clean, k = 10)
# convert to dataframe 
knn_clean <- as.data.frame(clean_kNN_imputed_30$data)
# display the sum of all NA values to check that they have all been imputed 
sum(is.na(knn_clean))
# print the first 6 rows to check the data has been imputed
head(knn_clean)
# Save the data frame as a text file
write.table(knn_clean, "knn_clean.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# missforest imputation
# run on matrix and not transposed 
mf_imputed_30 <- missForest(mat_countdata30_clean)
# convert the data with imputations to dataframe
mf <- as.data.frame(mf_imputed_30$ximp)
# check all NAs have been removed 
sum(is.na(mf))
# check it is formatted correctly
head(mf)
# save as a spreadsheet
write.table(mf, paste0(my_path, "mf_countdata30_imputed.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)

# hotdeck imputation
# hot-deck imputation function
hot_deck_impute <- function(data) {
  imputed_data <- data
  for (j in 1:ncol(data)) {
    missing <- is.na(data[, j])
    imputed_data[missing, j] <- mean(data[!missing, j], na.rm = TRUE)
  }
  
  return(imputed_data)
}
# run the hot deck imputation on the 30%
hot_deck <- hot_deck_impute(countdata30)
# check there is no NA values left
sum(is.na(hot_deck))
# save as a spreadsheet
write.table(hot_deck, paste0(my_path, "hotdeck_30.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
# inspect first 6 rows
head(hot_deck)

# DESeq2 analysis 
## Install BiocManager if not yet
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

# Install missing dependencies manually
install.packages("methods")
install.packages("tools")

# Ensure BiocManager is up-to-date
install.packages("BiocManager")

# Load BiocManager library
library(BiocManager)

# Remove the existing DESeq2 package if it exists
if ("DESeq2" %in% installed.packages()) {
  remove.packages("DESeq2")
}

# Use BiocManager to force re-install DESeq2 with dependencies
BiocManager::install("DESeq2", force = TRUE, dependencies = TRUE)

# Load DESeq2 library
library(DESeq2)
library(ggplot2) # To plot data
library(dplyr)
library(tibble)
library(gprofiler2)

BiocManager::install('EnhancedVolcano') # package for pretty volcano plots
library(EnhancedVolcano)


# Read in the count data and metadata
knn_30_countdata <- read.table(paste0(my_path, "knn_clean.txt"), header = TRUE, row.names = 1)


# Subset countdata to match coldata
knn_30_countdata <- knn_30_countdata[, colnames(knn_30_countdata) %in% coldata$RNAseq_Run_NCBI]

# Ensure the columns of countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(knn_30_countdata), coldata$RNAseq_Run_NCBI), ]

# Check if the columns match after subsetting
if (!all(colnames(knn_30_countdata) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}

# Remove rows with NA values
knn_30_countdata <- na.omit(knn_30_countdata)

# Round the count data to ensure all values are integers
knn_30_countdata <- round(knn_30_countdata)

# Create a DESeq object
knn_30_dds <- DESeqDataSetFromMatrix(countData = knn_30_countdata,
                                     colData = coldata,
                                     design = ~ Sex + group)

# Proceed with the differential expression analysis
knn_30_dds <- DESeq(knn_30_dds)

# Results
knn_30_res <- results(knn_30_dds)
knn_30_results_df <- as.data.frame(knn_30_res)
knn_30_results_df$Ensembl_id <- row.names(knn_30_results_df)
knn_30_results_df <- knn_30_results_df[order(knn_30_results_df$padj), ]

# Sanity check: Number of genes before the merge
knn_30_num_genes_before_merge <- nrow(knn_30_results_df)
cat("Number of genes before merge:", knn_30_num_genes_before_merge, "\n")

# Convert Ensembl names to gene names (if they are known)
knn_30_results_genes <- gconvert(row.names(knn_30_res), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)

# Ensure unique entries in results_genes
knn_30_results_genes <- knn_30_results_genes[!duplicated(knn_30_results_genes$input), ]

# Add the gene names
knn_30_results_df <- merge(knn_30_results_df,
                           knn_30_results_genes[, c("input", "target", "name", "description")],
                           by.x = "Ensembl_id", by.y = "input", all.x = TRUE)

# Sanity check: Number of genes after the merge
knn_30_num_genes_after_merge <- nrow(knn_30_results_df)
cat("Number of genes after merge:", knn_30_num_genes_after_merge, "\n")

# Ensure the 'Name' column is properly populated
knn_30_results_df$Name <- ifelse(is.na(knn_30_results_df$name), knn_30_results_df$Ensembl_id, knn_30_results_df$name)

# Check for missing values in the 'Name' column
knn_30_num_na_names <- sum(is.na(knn_30_results_df$Name))
cat("Number of missing names after merge:", knn_30_num_na_names, "\n")

# Subset the results to keep only significant genes
knn_30_significant_genes <- knn_30_results_df[knn_30_results_df$padj < 0.05 & !is.na(knn_30_results_df$padj), ]

# Ensure significant_genes is not empty before plotting
if (nrow(knn_30_significant_genes) > 0) {
  # Volcano plot
  EnhancedVolcano(knn_30_significant_genes,
                  lab = knn_30_significant_genes$Name,
                  x = 'log2FoldChange',
                  y = 'padj',
                  drawConnectors = TRUE)
} else {
  message("No significant genes found.")
}

# for missforest 
# Read in the count data and metadata
mf_countdata <- read.table(paste0(my_path, "mf_countdata30_imputed.txt"), header = TRUE, row.names = 1)

# Subset countdata to match coldata
mf_countdata <- mf_countdata[, colnames(mf_countdata) %in% coldata$RNAseq_Run_NCBI]

# Ensure the columns of countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(mf_countdata), coldata$RNAseq_Run_NCBI), ]

# Check if the columns match after subsetting
if (!all(colnames(mf_countdata) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}

# Remove rows with NA values
mf_countdata <- na.omit(mf_countdata)

# Round the count data to ensure all values are integers
mf_countdata <- round(mf_countdata)

# Create a DESeq object
mf_dds <- DESeqDataSetFromMatrix(countData = mf_countdata,
                                 colData = coldata,
                                 design = ~ Sex + group)

# Proceed with the differential expression analysis
mf_dds <- DESeq(mf_dds)

# Results
mf_res <- results(mf_dds)
mf_results_df <- as.data.frame(mf_res)
mf_results_df$Ensembl_id <- row.names(mf_results_df)
mf_results_df <- mf_results_df[order(mf_results_df$padj), ]

# Sanity check: Number of genes before the merge
mf_num_genes_before_merge <- nrow(mf_results_df)
cat("Number of genes before merge:", mf_num_genes_before_merge, "\n")

# Convert Ensembl names to gene names (if they are known)
mf_results_genes <- gconvert(row.names(mf_res), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)

# Ensure unique entries in results_genes
mf_results_genes <- mf_results_genes[!duplicated(mf_results_genes$input), ]

# Add the gene names
mf_results_df <- merge(mf_results_df,
                       mf_results_genes[, c("input", "target", "name", "description")],
                       by.x = "Ensembl_id", by.y = "input", all.x = TRUE)

# Sanity check: Number of genes after the merge
mf_num_genes_after_merge <- nrow(mf_results_df)
cat("Number of genes after merge:", mf_num_genes_after_merge, "\n")

# Ensure the 'Name' column is properly populated
mf_results_df$Name <- ifelse(is.na(mf_results_df$name), mf_results_df$Ensembl_id, mf_results_df$name)

# Check for missing values in the 'Name' column
mf_num_na_names <- sum(is.na(mf_results_df$Name))
cat("Number of missing names after merge:", mf_num_na_names, "\n")

# Subset the results to keep only significant genes
mf_significant_genes <- mf_results_df[mf_results_df$padj < 0.05 & !is.na(mf_results_df$padj), ]

# Ensure significant_genes is not empty before plotting
if (nrow(mf_significant_genes) > 0) {
  # Volcano plot
  EnhancedVolcano(mf_significant_genes,
                  lab = mf_significant_genes$Name,
                  x = 'log2FoldChange',
                  y = 'padj',
                  drawConnectors = TRUE)
} else {
  message("No significant genes found.")
}

# hotdeck DESeq2
# Read in the count data and metadata
hotdeck_countdata <- read.table(paste0(my_path, "hotdeck_30.txt"), header = TRUE, row.names = 1)

# Subset countdata to match coldata
hotdeck_countdata <- hotdeck_countdata[, colnames(hotdeck_countdata) %in% coldata$RNAseq_Run_NCBI]

# Ensure the columns of countdata correspond to the rows of coldata
coldata <- coldata[match(colnames(hotdeck_countdata), coldata$RNAseq_Run_NCBI), ]

# Check if the columns match after subsetting
if (!all(colnames(hotdeck_countdata) == coldata$RNAseq_Run_NCBI)) {
  stop("Count data columns do not match sample information rows.")
}

# Remove rows with NA values
hotdeck_countdata <- na.omit(hotdeck_countdata)

# Round the count data to ensure all values are integers
hotdeck_countdata <- round(hotdeck_countdata)

# Create a DESeq object
h_dds <- DESeqDataSetFromMatrix(countData = hotdeck_countdata,
                                colData = coldata,
                                design = ~ Sex + group)

# Proceed with the differential expression analysis
h_dds <- DESeq(h_dds)

# Results
h_res <- results(h_dds)
h_results_df <- as.data.frame(h_res)
h_results_df$Ensembl_id <- row.names(h_results_df)
h_results_df <- h_results_df[order(h_results_df$padj), ]

# Sanity check: Number of genes before the merge
h_num_genes_before_merge <- nrow(h_results_df)
cat("Number of genes before merge:", h_num_genes_before_merge, "\n")

# Convert Ensembl names to gene names (if they are known)
h_results_genes <- gconvert(row.names(h_res), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)

# Ensure unique entries in results_genes
h_results_genes <- h_results_genes[!duplicated(h_results_genes$input), ]

# Add the gene names
h_results_df <- merge(h_results_df,
                      h_results_genes[, c("input", "target", "name", "description")],
                      by.x = "Ensembl_id", by.y = "input", all.x = TRUE)

# Sanity check: Number of genes after the merge
h_num_genes_after_merge <- nrow(h_results_df)
cat("Number of genes after merge:", h_num_genes_after_merge, "\n")

# Ensure the 'Name' column is properly populated
h_results_df$Name <- ifelse(is.na(h_results_df$name), h_results_df$Ensembl_id, h_results_df$name)

# Check for missing values in the 'Name' column
h_num_na_names <- sum(is.na(h_results_df$Name))
cat("Number of missing names after merge:", h_num_na_names, "\n")

# Subset the results to keep only significant genes
h_significant_genes <- h_results_df[h_results_df$padj < 0.05 & !is.na(h_results_df$padj), ]

# Ensure significant_genes is not empty before plotting
if (nrow(h_significant_genes) > 0) {
  # Volcano plot
  EnhancedVolcano(h_significant_genes,
                  lab = h_significant_genes$Name,
                  x = 'log2FoldChange',
                  y = 'padj')
} else {
  message("No significant genes found.")
}

# EVALUATION OF METHODS
# knn
knn_imputation_evaluation(paste0(my_path, "knn_clean.txt"))

# missforest 
missforest_evaluation_function(mf_countdata)


# hotdeck 
hotdeck_imputation_evaluation <- function(file_path = paste0(my_path, "hotdeck_30.txt"), na_fraction = 0.3, seed1 = 123, seed2 = 456) {
  # Load the dataset
  countdata30_imputedHotDeck <- read.table(file_path, header = TRUE, row.names = 1)
  
  # Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed1) # For reproducibility
  total_samples <- nrow(countdata30_imputedHotDeck)
  train_size <- floor(0.7 * total_samples)
  
  train_indices <- sample(seq_len(total_samples), size = train_size)
  train_data <- countdata30_imputedHotDeck[train_indices, ]
  test_data <- countdata30_imputedHotDeck[-train_indices, ]
  
  # Introduce artificial NAs randomly in the testing set
  set.seed(seed2) # For reproducibility
  test_data_with_na <- test_data
  num_values <- prod(dim(test_data_with_na))
  num_na <- round(num_values * na_fraction)
  na_indices <- arrayInd(sample(num_values, num_na), dim(test_data_with_na))
  test_data_with_na[na_indices] <- NA
  
  # Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.")
  }
  
  # Re-impute the missing values in the testing set using hot deck
  test_data_reimputed <- hot_deck_impute(test_data_with_na)$data
  
  # Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.")
  }
  
  # Extract true and imputed values for comparison
  true_values <- test_data[na_indices]
  imputed_values <- test_data_reimputed[na_indices]
  
  # Check if there are differences between true and imputed values
  if (all(true_values == imputed_values)) {
    cat("No differences between true and imputed values.\n")
    mae_value <- 0
    rmse_value <- 0
  } else {
    # Calculate MAE and RMSE between the imputed values and the true values in the test set
    mae_value <- mae(true_values, imputed_values)
    rmse_value <- rmse(true_values, imputed_values)
  }
  
  cat("Mean Absolute Error (MAE):", mae_value, "\n")
  cat("Root Mean Squared Error (RMSE):", rmse_value, "\n")
  
  # Save the datasets
  write.table(train_data, file = paste0(my_path, "countdata30_train.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_with_na, file = paste0(my_path, "countdata30_test_with_na.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  write.table(test_data_reimputed, file = paste0(my_path, "countdata30_test_reimputed.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Plot the distribution of true values and imputed values
  true_vs_imputed <- data.frame(True = true_values, Imputed = imputed_values)
  
  # Plot a scatter plot and a line y=x for reference
  p1 <- ggplot(true_vs_imputed, aes(x = True, y = Imputed)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    labs(title = "True vs Imputed Values", x = "True Values", y = "Imputed Values") +
    theme_minimal()
  
  # Plot the error distribution
  errors <- data.frame(Error = true_values - imputed_values)
  
  p2 <- ggplot(errors, aes(x = Error)) +
    geom_histogram(binwidth = 10, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Error Distribution (True - Imputed)", x = "Error", y = "Frequency") +
    theme_minimal()
  
  print(p1)
  print(p2)
}

# use function on the hotdeck imputed dataframe
hotdeck_imputation_evaluation(paste0(my_path, "hotdeck_30.txt"))
