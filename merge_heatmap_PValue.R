# Load necessary libraries
setwd("D:/rnaseq_exp/exp9/pvalue")
# Load the CSV file (no 'Description' column to handle)
cluster_matrix_hm = read.csv("z_score_PValue_all.csv", stringsAsFactors = F, row.names = 1)

# Create merged data frame by averaging samples 32 and 33, and samples 34 and 35
cluster_matrix_merged = data.frame(cbind(  
  rowMeans(cluster_matrix_hm[,1:3]),  # Merge sample 32 and 33
  rowMeans(cluster_matrix_hm[,2:4])   # Merge sample 34 and 35
), stringsAsFactors = F)

# Rename columns accordingly
colnames(cluster_matrix_merged) = c("Sample32_33", "Sample34_35")

# Normalize and scale the data matrix
matrix_norm = t(scale(t(cluster_matrix_merged), center = T, scale = T))
matrix_norm[is.na(matrix_norm)] <- 0

# Define new column names for heatmap
cn = c("Sample32_33", "Sample34_35")

# Create heatmap
library(ComplexHeatmap)
hm = Heatmap(matrix_norm, 
             cluster_rows = T, name = " ", 
             row_names_gp = gpar(fontsize = 12), show_column_dend = F,
             show_row_names = T,
             column_names_gp = gpar(fontsize=12, fontface="italic"),
             bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 45, just = "right", gp=gpar(fontface="italic")),annotation_height = max_text_width(cn)),
             show_row_dend = F,
             show_column_names = F,
             left_annotation = rowAnnotation(pt = anno_empty(border = F, width = unit(5,"mm")))
)

# Draw the heatmap
draw(hm, padding = unit(c(2, 0, 10, 2), "mm"))

