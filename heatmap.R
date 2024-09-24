# Load necessary libraries
library(ComplexHeatmap)
library(circlize)  # For colorRamp

# Set working directory
setwd("D:/rnaseq_exp/exp10")

# Load your data from CSV
data <- read.csv("z_score_log2.csv", row.names = 1)

# Cap Z-scores between -2 and +2 (if not already done)
data[data > 2] <- 2
data[data < -2] <- -2

# Perform clustering
row_dist <- dist(data)
row_clustering <- hclust(row_dist, method = "complete")
row_clusters <- cutree(row_clustering, k = 10)  # Adjust 'k' as needed

col_dist <- dist(t(data))
col_clustering <- hclust(col_dist, method = "complete")
col_clusters <- cutree(col_clustering, k = 2)  # Adjust 'k' as needed

# Create annotation data frame for rows
row_annotation <- data.frame(Cluster = factor(row_clusters))



# Convert annotations to HeatmapAnnotation objects
row_annotation_ha <- rowAnnotation(df = row_annotation, 
                                   show_legend = TRUE,  # Show legend for row clusters
                                   annotation_legend_param = list(title = "Clusters"))
# Save the result to a new CSV file
#write.csv(row_annotation, "cluster_results_exp-10.csv", row.names = TRUE)

# Define color palette for heatmap
color_palette <- colorRampPalette(c("red", "white", "blue"))(5000)

# Generate the heatmap
Heatmap(data,
        name = "Z-score",
        col = color_palette,
        cluster_rows = row_clustering,
        cluster_columns = col_clustering,
        show_row_names = FALSE,
        show_column_names = TRUE,  # Hide column names
        right_annotation = row_annotation_ha,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        show_heatmap_legend = TRUE,  # Hide the color palette legend
        heatmap_legend_param = list(title = "Z-score"))  # Set the title if needed

