# Set the working directory (change this to your folder containing CSV files)
setwd("D:/rnaseq_exp/GO David/exp8/analysis/csv")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Get the list of all CSV files in the directory
file_list <- list.files(pattern = "\\.csv$")

# Loop through each CSV file
for (file in file_list) {
  
  # Load the GO results for the current file
  go_results <- read.csv(file)
  
  # Filter terms with PValue <= 0.05
  filtered_terms <- go_results %>% filter(PValue <= 0.05)
  
  
  # Filter terms with PValue <= 0.05 and take the top 30 based on Count
  #filtered_terms <- go_results %>%
   # filter(PValue <= 0.05) %>%
    #top_n(30, Count)
  
  
  # Extract the base name (without extension) for dynamic naming
  base_name <- tools::file_path_sans_ext(file)
  
  # Create a barplot for the filtered GO terms by gene count with color gradient based on PValue
  barplot <- ggplot(filtered_terms, aes(x = reorder(Term, -Count), y = Count, fill = PValue)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue") +  # Color by p-value
    coord_flip() +  # Horizontal bars
    theme_minimal() +  # Basic minimal theme
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background for graph area
      plot.background = element_rect(fill = "white", color = NA),   # White background for complete image
      panel.grid.major = element_line(color = "grey90"),            # Light grid lines
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis font size
      axis.text.y = element_text(size = 14),                        # Increase y-axis font size
      axis.title.x = element_text(size = 16, face = "bold"),        # Increase x-axis title font size
      axis.title.y = element_text(size = 16, face = "bold"),        # Increase y-axis title font size
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Increase plot title size
    ) +
    labs(title = paste("GO Terms (PValue <= 0.05) by Gene Count -", base_name), 
         x = "GO Terms", y = "Gene Count", fill = "PValue")
  
  # Save the barplot with dynamic file name
  ggsave(paste0(base_name, "_Barplot_PValue_0.05.png"), plot = barplot, width = 20, height = 12, bg = "white")
  
  # Create a dotplot for the filtered GO terms with a white background and p-value color gradient
  dotplot <- ggplot(filtered_terms, aes(x = Fold.Enrichment, y = reorder(Term, Fold.Enrichment))) +
    geom_point(aes(size = Count, color = PValue)) +
    scale_color_gradient(low = "red", high = "blue") +  # Color by p-value
    theme_minimal() +  # Basic minimal theme
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background for graph area
      plot.background = element_rect(fill = "white", color = NA),   # White background for complete image
      panel.grid.major = element_line(color = "grey90"),            # Light grid lines
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis font size
      axis.text.y = element_text(size = 14),                        # Increase y-axis font size
      axis.title.x = element_text(size = 16, face = "bold"),        # Increase x-axis title font size
      axis.title.y = element_text(size = 16, face = "bold"),        # Increase y-axis title font size
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # Increase plot title size
    ) +
    labs(title = paste("GO Terms Dotplot: Fold Enrichment vs Gene Count -", base_name), 
         x = "Fold Enrichment", y = "GO Terms", color = "PValue", size = "Gene Count")
  
  # Save the dotplot with dynamic file name
  ggsave(paste0(base_name, "_Dotplot_PValue_0.05.png"), plot = dotplot, width = 20, height = 12, bg = "white")
}

