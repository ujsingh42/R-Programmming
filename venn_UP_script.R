# Load necessary libraries
library(VennDiagram)
library(grid)

# Set working directory
setwd("D:/rnaseq_exp/ven digram")

# Read the CSV file and clean empty values (empty strings and spaces)
df <- read.csv("UP.csv", stringsAsFactors = FALSE)

# Replace empty strings and spaces in the data frame with NA
df[df == "" | df == " "] <- NA

# Create sets for siMECOM KD and JIB04 treatment
set_exp9 <- df$exp9
set_exp10 <- df$exp10

# Remove NA values from sets to avoid errors in the Venn diagram
set_exp9 <- na.omit(set_exp9)
set_exp10 <- na.omit(set_exp10)

# Identify and count duplicates in set_exp9
duplicates_exp9 <- set_exp9[duplicated(set_exp9)]
cat("Number of duplicate values in set_exp9:", length(duplicates_exp9), "\n")

# Remove duplicates from set_exp9
set_exp9 <- unique(set_exp9)
cat("Number of unique values in set_exp9 after removing duplicates:", length(set_exp9), "\n")

# Identify the intersection (overlap) and differences between sets
overlap <- intersect(set_exp9, set_exp10)
only_exp9 <- setdiff(set_exp9, set_exp10)
only_exp10 <- setdiff(set_exp10, set_exp9)

# Calculate the total count for siMECOM set
total_siMECOM <- length(only_exp9) + length(overlap)  # Unique + Overlap
cat("Correct total count for siMECOM set:", total_siMECOM, "\n")

# Check for any empty strings or blanks in only_exp9 and replace them with NA
only_exp9[only_exp9 == "" | only_exp9 == " "] <- NA

# Create a data frame to store the information and replace empty values with NA
max_len <- max(length(only_exp9), length(only_exp10), length(overlap))

# Adjust the lengths of vectors by padding with NAs where needed
venn_data <- data.frame(
  "siMECOM_KD_only" = c(only_exp9, rep(NA, max_len - length(only_exp9))),
  "JIB04_upreg_only" = c(only_exp10, rep(NA, max_len - length(only_exp10))),
  "Overlap" = c(overlap, rep(NA, max_len - length(overlap)))
)

# Save the results in a CSV file without blank rows, replacing them with NA
write.csv(venn_data, file = "venn_diagram_down_results.csv", row.names = FALSE, na = "NA")

# Create a Venn diagram with adjusted text size and increased background margins
venn.plot <- venn.diagram(
  x = list("OVSAHO siMECOM" = set_exp9, "OVSAHO JIB04-upreg" = set_exp10),
  filename = NULL,                     # We will draw it to an object first
  fill = c("#FF9999", "#9999FF"),       # Custom pastel colors
  alpha = c(0.7, 0.7),                  # Transparency to see overlaps better
  lty = "solid",                        # Line type: solid
  lwd = 2,                              # Line width: thicker lines
  col = c("#FF6666", "#6666FF"),        # Border colors
  cex = 1.8,                            # Reduced text size for counts (inside the circles)
  fontface = "bold",                    # Bold text for counts
  cat.cex = 3,                          # Increased size of set labels
  cat.fontface = "bold",                # Bold font for set labels
  cat.default.pos = "outer",            # Position labels outside the circles
  cat.pos = c(10, 0.50),                   # Adjust position of labels
  cat.dist = c(0.05, 0.05),             # Distance from the circles
  cat.col = c("#FF6666", "#6666FF"),    # Label colors
  main = "OVSAHO JIB04-upreg vs OVSAHO siMECOM",  # Add a title
  main.cex = 4,                         # Title font size
  main.fontface = "bold",               # Title in bold
  sub = "Gene Set Overlap",             # Add a subtitle
  sub.cex = 4,                           # Subtitle font size
  scaled = FALSE                        # Set to FALSE to make both circles the same size
  )

# Save the Venn diagram to a PNG file with larger background margins
png("OVSAHO JIB04-upreg vs OVSAHO siMECOM.png", width = 1600, height = 1600, res = 100)  # Increased canvas size
par(mar = c(6, 6, 6, 6))  # Adjust margins: bottom, left, top, right
grid.draw(venn.plot)
dev.off()

# Save the Venn diagram to a PDF file with larger background margins
pdf("OVSAHO JIB04-upreg vs OVSAHO siMECOM.pdf", width = 16, height = 16)  # Larger canvas
par(mar = c(6, 6, 6, 6))  # Adjust margins: bottom, left, top, right
grid.draw(venn.plot)
dev.off()

# Draw the Venn diagram in R's graphical window
grid.draw(venn.plot)

