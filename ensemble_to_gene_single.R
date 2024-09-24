# Load necessary libraries
library(biomaRt)
library(readr)

# Load necessary libraries
setwd("D:/rnaseq_exp/GO David/exp10")

# Define the input and output file paths
input_file <- "pvalue_exp10.csv"  # Ensure the file extension is included
output_file <- "pvalue_exp10_reuslt_Genes.csv"

# Load Ensembl IDs from the CSV file
ensembl_data <- read.csv(input_file)

# Print the first few rows to verify the data
print(head(ensembl_data))

# Check if the column with Ensembl IDs is named "ensembl_id"
# If not, replace "ensembl_id" with the correct column name in the following line
ensembl_ids <- ensembl_data$ensembl_id

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query to convert Ensembl IDs to gene symbols
result <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Merge the result with the original data to keep the Ensembl IDs
final_result <- merge(
  ensembl_data,
  result,
  by.x = "ensembl_id",  # Column name in original data
  by.y = "ensembl_gene_id",  # Column name in Ensembl results
  all.x = TRUE
)

# Save the result to a new CSV file
write_csv(final_result, output_file)

# Print the first few rows of the final result to verify
print(head(final_result))


