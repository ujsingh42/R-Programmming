# Set working directory
setwd("D:/rnaseq_exp/GO David/exp9/analysis/csv")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load DAVID GO results
david_go_results <- read.csv("BP_down.csv")

# Check the structure of the data
str(david_go_results)

# Select the top 30 terms based on Count
top_terms <- david_go_results %>% top_n(30, Count)

# Create a barplot for the top 30 GO terms by gene count with color gradient based on PValue
barplot <- ggplot(top_terms, aes(x = reorder(Term, -Count), y = Count, fill = PValue)) +
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
  labs(title = "Top 30 GO Terms by Gene Count (BP Downregulated)", 
       x = "GO Terms", y = "Gene Count", fill = "PValue")

# Save the barplot with white background
ggsave("GO_BP_Downregulated_Barplot.png", plot = barplot, width = 20, height = 12, bg = "white")

# Create a dotplot for the top 30 GO terms with a white background and p-value color gradient
dotplot <- ggplot(top_terms, aes(x = Fold.Enrichment, y = reorder(Term, Fold.Enrichment))) +
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
  labs(title = "GO Terms Dotplot: Fold Enrichment vs Gene Count (BP Downregulated)", 
       x = "Fold Enrichment", y = "GO Terms", color = "PValue", size = "Gene Count")

# Save the dotplot with white background
ggsave("GO_BP_Downregulated_Dotplot.png", plot = dotplot, width = 20, height = 12, bg = "white")

