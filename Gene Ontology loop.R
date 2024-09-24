# Load necessary libraries
setwd("D:/rnaseq_exp/exp8/GO")
library(clusterProfiler)
library(org.Hs.eg.db)  # For human gene annotations
library(DOSE)
library(enrichplot)
library(ggplot2)

# Import edgeR results CSV
edger_results <- read.csv("results_exp8_gene_symbols.csv", row.names = 1)

# Filter genes based on p-value < 0.05
filtered_DEGs <- edger_results[edger_results$PValue < 0.05, ]

# Split into upregulated and downregulated genes based on logFC
upregulated_genes <- rownames(filtered_DEGs[filtered_DEGs$logFC > 0, ])
downregulated_genes <- rownames(filtered_DEGs[filtered_DEGs$logFC < 0, ])

# Convert upregulated gene identifiers to Entrez IDs
gene_list_up <- bitr(upregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez_ids_up <- gene_list_up$ENTREZID

# Convert downregulated gene identifiers to Entrez IDs
gene_list_down <- bitr(downregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez_ids_down <- gene_list_down$ENTREZID

# Perform GO enrichment analysis for upregulated genes
ego_up_BP <- enrichGO(gene = gene_entrez_ids_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05)
ego_up_CC <- enrichGO(gene = gene_entrez_ids_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05)
ego_up_MF <- enrichGO(gene = gene_entrez_ids_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05)

# Perform GO enrichment analysis for downregulated genes
ego_down_BP <- enrichGO(gene = gene_entrez_ids_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05)
ego_down_CC <- enrichGO(gene = gene_entrez_ids_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05)
ego_down_MF <- enrichGO(gene = gene_entrez_ids_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05)

# Plot and save results for each GO term category (BP, CC, MF)

# Upregulated genes - BP, CC, MF
barplot(ego_up_BP, showCategory = 40, title = "GO Enrichment for BP Upregulated Genes")
ggsave("BP_UP_GO_Plot.png", width = 8, height = 12)

barplot(ego_up_CC, showCategory = 40, title = "GO Enrichment for CC Upregulated Genes")
ggsave("CC_UP_GO_Plot.png", width = 8, height = 12)

barplot(ego_up_MF, showCategory = 40, title = "GO Enrichment for MF Upregulated Genes")
ggsave("MF_UP_GO_Plot.png", width = 8, height = 12)

# Downregulated genes - BP, CC, MF
barplot(ego_down_BP, showCategory = 40, title = "GO Enrichment for BP Downregulated Genes")
ggsave("BP_DOWN_GO_Plot.png", width = 8, height = 12)

barplot(ego_down_CC, showCategory = 40, title = "GO Enrichment for CC Downregulated Genes")
ggsave("CC_DOWN_GO_Plot.png", width = 8, height = 12)

barplot(ego_down_MF, showCategory = 40, title = "GO Enrichment for MF Downregulated Genes")
ggsave("MF_DOWN_GO_Plot.png", width = 8, height = 12)

# Save the results in CSV format for each category
write.csv(as.data.frame(ego_up_BP), file = "GO_BP_Upregulated_Results.csv")
write.csv(as.data.frame(ego_up_CC), file = "GO_CC_Upregulated_Results.csv")
write.csv(as.data.frame(ego_up_MF), file = "GO_MF_Upregulated_Results.csv")

write.csv(as.data.frame(ego_down_BP), file = "GO_BP_Downregulated_Results.csv")
write.csv(as.data.frame(ego_down_CC), file = "GO_CC_Downregulated_Results.csv")
write.csv(as.data.frame(ego_down_MF), file = "GO_MF_Downregulated_Results.csv")
