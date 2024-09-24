# Load necessary libraries
setwd("D:/rnaseq_exp/exp9/pvalue")
#install.packages("openxlsx")
#install.packages("BiocManager")
#BiocManager::install("ggfortify")
library(edgeR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggfortify)
library(ComplexHeatmap)
library(circlize)
library(readr)

# Load count data
counts <- read.csv("count.csv", row.names = 1)

# Load metadata
metadata <- read.csv("metadata.csv")

# Create a DGEList object
group <- factor(metadata$condition)
dge <- DGEList(counts = counts, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize data
dge <- calcNormFactors(dge)

# Compute log-transformed CPM values
logCPM <- cpm(dge, log = TRUE)

# PCA plot
autoplot(prcomp(t(logCPM)), data = metadata, colour = 'condition', shape = 17) + # shape = 16 corresponds to filled circles
  theme_minimal() +
  labs(title = "PCA Plot ")

# Create design matrix
design <- model.matrix(~group)

# Estimate dispersion
dge <- estimateDisp(dge, design)
plotBCV(dge)

# Fit the model and perform statistical tests
fit <- glmQLFit(dge, design)
plotQLDisp(fit)
qlf <- glmQLFTest(fit, coef = 2) # Change coef based on comparison
plotMD(qlf)


# Get results
results <- topTags(qlf, n = Inf)

# Save results
write.csv(results$table, "edgeR_exp9_reuslt.csv")

# MA Plot
plotMD(qlf)
abline(h = c(-1, 1), col = "red")


########### Volcano Plot FDR < 0.05

# Assuming `results$table` is your data frame containing the results
results_table <- results$table

# Add significance column
results_table$Significance <- ifelse(results_table$FDR < 0.05, "Significant", "Not Significant")

# Calculate -log10(FDR) for the y-axis
results_table$neg_log10_FDR <- -log10(results_table$FDR)

# Classify genes as Upregulated, Downregulated, or Not Significant
results_table$Regulation <- with(results_table, 
                                 ifelse(FDR < 0.05 & logFC > 1, "Upregulated",
                                        ifelse(FDR < 0.05 & logFC < -1, "Downregulated", "Not Significant")))

# Create the volcano plot with brighter colors
ggplot(results_table, aes(x = logFC, y = neg_log10_FDR, color = Regulation)) +
  geom_point(alpha = 0.8, size = 2) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log Fold Change (logFC)", y = "-log10(FDR)") +
  scale_color_manual(values = c("Upregulated" = "#FF0000", "Downregulated" = "#0000FF", "Not Significant" = "#000000")) +
  theme(legend.position = "top") # Optional: place legend at the top

################ Volcano Plot PValue < 0.05

# Assuming `results$table` is your data frame containing the results
results_table <- results$table

# Add significance column
results_table$Significance <- ifelse(results_table$P < 0.05, "Significant", "Not Significant")

# Calculate -log10(FDR) for the y-axis
results_table$neg_log10_FDR <- -log10(results_table$FDR)

# Classify genes as Upregulated, Downregulated, or Not Significant
results_table$Regulation <- with(results_table, 
                                 ifelse(PValue < 0.05 & logFC > 1, "Upregulated",
                                        ifelse(PValue < 0.05 & logFC < -1, "Downregulated", "Not Significant")))

# Filter data to include only significant results
filtered_results <- subset(results_table, PValue < 0.05)

# Create the volcano plot
ggplot(filtered_results, aes(x = logFC, y = neg_log10_FDR, color = Regulation)) +
  geom_point(alpha = 0.8, size = 2) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log Fold Change (logFC)", y = "-log10(FDR)") +
  scale_color_manual(values = c("Upregulated" = "#FF0000", "Downregulated" = "#0000FF", "Not Significant" = "#000000")) +
  theme(legend.position = "top") # Optional: place legend at the top


############## Z score from CPM PValue < 0.05

# Filter DEGs based on PVAlue < 0.05
filtered_DEGs <- results$table[results$table$PValue < 0.05, ]


# Extract gene names of the filtered DEGs
filtered_genes <- rownames(filtered_DEGs)

# Retrieve counts for the filtered DEGs
filtered_counts <- cpm(dge)[filtered_genes, ]

# Calculate Z-scores for each gene (row)
z_scores_all <- t(scale(t(filtered_counts)))

# Convert matrix to data frame and include gene names
z_scores_df <- as.data.frame(z_scores_all)
z_scores_df$Gene <- rownames(z_scores_all)  # Add gene names as a column

# Reorder columns to have 'Gene' first
z_scores_df <- z_scores_df[, c("Gene", setdiff(names(z_scores_df), "Gene"))]

# Write the data frame to a CSV file with gene names
write.csv(z_scores_df, file = "z_score_PValue_all.csv", row.names = FALSE)

# Cap Z-scores between -2 and +2
z_scores_all[z_scores_all > 2] <- 2
z_scores_all[z_scores_all < -2] <- -2

# Convert matrix to data frame and include gene names
z_scores_capped_df <- as.data.frame(z_scores_all)
z_scores_capped_df$Gene <- rownames(z_scores_all)  # Add gene names as a column

# Reorder columns to have 'Gene' first
z_scores_capped_df <- z_scores_capped_df[, c("Gene", setdiff(names(z_scores_capped_df), "Gene"))]

# Write the capped Z-scores data frame to a CSV file
write.csv(z_scores_capped_df, file = "z_score_PValue_cap_to_2.csv", row.names = FALSE)

# Load your data from CSV
data <- read.csv("z_score_PValue_all.csv", row.names = 1)

# Cap Z-scores between -2 and +2 
data[data > 2] <- 2
data[data < -2] <- -2

# Convert data to a matrix
data_matrix <- as.matrix(data)

# Perform clustering on rows
row_dist <- dist(data_matrix)
row_clustering <- hclust(row_dist, method = "complete")
row_clusters <- cutree(row_clustering, k = 10)  # Adjust 'k' as needed

# Perform clustering on columns
col_dist <- dist(t(data_matrix))
col_clustering <- hclust(col_dist, method = "complete")
col_clusters <- cutree(col_clustering, k = 2)  # Adjust 'k' as needed

# Map clusters to "Treated" and "Control" labels
col_labels <- ifelse(col_clusters == 2, "Treated", "Control")

# Create annotation data frame for rows
row_annotation <- data.frame(Cluster = factor(row_clusters))

# Convert annotations to HeatmapAnnotation objects for rows
row_annotation_ha <- rowAnnotation(df = row_annotation, 
                                   show_legend = TRUE,  # Show legend for row clusters
                                   annotation_legend_param = list(title = "Clusters"))

# Create annotation data frame for columns
col_annotation <- data.frame(Treatment = factor(col_labels))

# Convert annotations to HeatmapAnnotation objects for columns
col_annotation_ha <- HeatmapAnnotation(df = col_annotation, 
                                       show_legend = TRUE,  # Show legend for column clusters
                                       annotation_legend_param = list(title = "Treatment"))

# Save the row cluster result to a new CSV file
write.csv(row_annotation, "cluster_results_exp-10.csv", row.names = TRUE)

# Define color palette for heatmap
color_palette <- colorRampPalette(c("blue", "white", "red"))(5000)

# Generate the heatmap with original column order
Heatmap(data_matrix,  # Use the original data matrix without reversing
        name = "Z-score",
        col = color_palette,
        cluster_rows = row_clustering,
        cluster_columns = col_clustering,
        show_row_dend = F, show_column_dend = F,  # Hide dendrograms
        show_row_names = FALSE,
        show_column_names = F,  # Show column names
        right_annotation = row_annotation_ha,
        top_annotation = col_annotation_ha,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        show_heatmap_legend = TRUE,  # Show the color palette legend
        heatmap_legend_param = list(title = "Z-score"))  # Set the title if needed

