library(DESeq2)
library(dplyr)
library(ggplot2)

# Load the counts data (assuming tab-delimited format)
counts <- read.table("data/processed_data/count_files/all_samples_counts.txt", sep="\t", header=TRUE, row.names=1)

# Check the structure of the counts data (first few rows)
head(counts)

# Remove the non-sample columns (Chr, Start, End, Strand, Length) and keep only the counts
counts <- counts[, -c(1:5)]  # Remove the first 5 columns (Chr, Start, End, Strand, Length)

# Check the modified counts data (should now only contain sample columns)
head(counts)

# Rename the column names to just the sample identifiers (e.g., ERR2713020, ERR2713021, etc.)
colnames(counts) <- gsub(".*sorted_(ERR\\d+).*", "\\1", colnames(counts))

# Check the new column names
print(colnames(counts))

# Handle duplicated row names by making them unique, then set as rownames in the dataframe
uniq_name <- make.names(rownames(counts), unique = TRUE)
rownames(counts) <- uniq_name

# Check the first few rows after modification
head(counts)

# Create sample information (you can adjust the experimental design as needed)
# Assuming you have 5 samples, each with the experimental groups "ctrl" and "wt" (adjust this based on your data)
samples <- data.frame(
  row.names = colnames(counts),
  theSample = rep(c("ctrl", "wt"), length.out = length(colnames(counts)))  # Adjust as necessary
)

# Ensure that the number of samples matches the number of columns in the count data
if (length(colnames(counts)) != nrow(samples)) {
  stop("Mismatch between the number of columns in count data and the number of sample labels")
}

# Create DESeq2 dataset object from the count matrix and sample information
ds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~ theSample)

# Ensure that column names in the DESeqDataSet match those in the count matrix
colnames(ds) <- colnames(counts)

# Run DESeq2 analysis
ds <- DESeq(ds)

# Extract results from DESeq2 analysis
res <- results(ds)

# Display the results
print(res)

# Identify significantly differentially expressed genes (adjusted p-value < 0.05)
sig <- res[which(res$padj < 0.05), ]
print(sig)

# Identify genes with significant log2 fold change (> 2 or < -2)
sigLf <- res[which(res$log2FoldChange < -2 | res$log2FoldChange > 2), ]
print(sigLf)

# Define the file path and name for saving the MA plot
png("my_MA_plot.png")

# Generate the MA plot
plotMA(ds)

# Close the device to save the plot
dev.off()


