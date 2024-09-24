# Load necessary libraries
setwd("D:/rnaseq_exp/exp6")
# Load necessary libraries
library(DESeq2)
library(readr)

# Load count data
counts <- read.csv("count.csv", row.names = 1)


# Load experimental design data
colData <- read.csv("metadata.csv", row.names = 1)

# Trim leading and trailing spaces from row names and column names
rownames(colData) <- trimws(rownames(colData))
colnames(counts) <- trimws(colnames(counts))

# Check that the sample names match
all(rownames(colData) == colnames(counts))


print(rownames(colData))
print(colnames(counts))



# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

# Run the DESeq function
dds <- DESeq(dds)
resHTSeq <- results(dds)

# Get results for the contrast of interest
all_genes_results <- results(dds, contrast = c("condition", "treated", "control"))
library(dplyr)
# distribution of adjusted p-values
hist(all_genes_results$padj, col="lightblue", main = "Adjusted p-value distribution")


resLFC <- lfcShrink(dds = dds, 
                    res = all_genes_results,
                    type = "normal",
                    coef = "condition_treated_vs_control") # name or number of the coefficient (LFC) to shrink



resultsNames(dds)

# load the library if not done yet
library("EnhancedVolcano")

# The main function is named after the package
EnhancedVolcano(toptable = resLFC,              # We use the shrunken log2 fold change as noise associated with low count genes is removed 
                x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
                y = "padj",                     # Name of the column in resLFC that contains the p-value
                lab = rownames(resLFC)
)



EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC),
                xlim = c(-10, 10),
                ylim = c(0, 100),
                pCutoff = 1e-06,
                pointSize = 2.0,
                FCcutoff = 2, 
                title = "condition_treated_vs_normal \n (fold change cutoff = 2, p-value cutoff = 1e-06)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)


head(resHTSeq)
table(resHTSeq$padj < 0.01)
orderedRes <- resHTSeq[ order(resHTSeq$pvalue), ]

write.csv(as.data.frame(orderedRes), file="deseq2_result_exp10.csv")
normCounts <- counts(dds, normalized = TRUE)

head(normCounts)
#write.csv(as.data.frame(orderedRes), file="exp2_normCounts.DESeq2.csv")
plotDispEsts(dds)
hist(resHTSeq$pvalue, breaks=0:50/50, xlab="p value", main="Histogram of nominal p values")
plotMA(resHTSeq)

plot(resHTSeq$log2FoldChange, -log10(resHTSeq$pvalue), xlab="log2 Fold-change", ylab="-log P-value", pch=20, cex=0.5)
points(resHTSeq$log2FoldChange[ resHTSeq$padj<0.05 ], -log10(resHTSeq$pvalue[ resHTSeq$padj<0.05 ]), col="red", pch=20, cex=0.5)
abline(v=0, h=-log10(0.05), lty="dashed", col="grey")


vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

plotPCA(vsd)

dists <- dist(t(assay(vsd)))
# headmap of distances
heatmap(as.matrix(dists), main="Clustering of euclidean distances", scale="none")

#install.packages("gplots")

library(gplots)

diffgenes <- rownames(resHTSeq)[ which(resHTSeq$pvalue < 0.05) ]
diffcounts <- normCounts[ diffgenes, ]

# Save as a PDF
pdf("heatmap.pdf", width = 10, height = 8)
heatmap.2(diffcounts, 
          labRow = "", 
          trace = "none", 
          density.info = "none",
          scale = "row",
          distfun = function(x) as.dist(1 - cor(t(x))))
dev.off()  # Close the PDF device

# Save as a PNG
png("heatmap.png", width = 1000, height = 800)
heatmap.2(diffcounts, 
          labRow = "", 
          trace = "none", 
          density.info = "none",
          scale = "row",
          distfun = function(x) as.dist(1 - cor(t(x))))
dev.off()  # Close the PNG device


