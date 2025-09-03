# conda activate r_env
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(cowplot)
library(ggplot2)
library(dplyr) 
library(tidyr)
library(stringr)

setwd("P1_Supertype_DE")
de <- readRDS("Files/DE_results_Affected_vs_Unaffected.rds")
sig_de <- subset(de, adj.P.Val < 0.05)
sig_de_alpha <- sig_de[order(sig_de$genes), ]

#Save just sig genes
write.csv(sig_de_alpha, "Files/DE_sig_genes_sst.csv", row.names = FALSE)


# Go enrichment
ego <- enrichGO(gene          = sig_de$genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "BP",          # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)



p <- barplot(ego, showCategory = 15, title = "GO: Affected v Unaffected")
ggsave("Figures/GO_Affected_vs_Unaffected2.png", p, width = 14, height = 10, dpi = 300)

p <- barplot(ego, showCategory = 8)
ggsave("Figures/GO_Affected_vs_Unaffected3.png", p, width = 14, height = 10, dpi = 300)
ggsave("Figures/GO_Affected_vs_Unaffected3.svg", p, width = 14, height = 10, dpi = 300)



# Convert to dataframe
ego_df <- as.data.frame(ego)


terms <- c( "regulation of synapse structure or activity",
  "regulation of synapse organization",
  "regulation of membrane potential")

selected_terms <- ego_df[ego_df$Description %in% terms, ]

# Expand into a tidy dataframe of term-gene pairs
term_gene_list <- lapply(1:nrow(selected_terms), function(i) {
  genes <- unlist(strsplit(selected_terms$geneID[i], "/"))
  data.frame(Term = selected_terms$Description[i], genes = genes, stringsAsFactors = FALSE)
}) %>% bind_rows()

# Merge with DE results (logFC + adj.P.Val)
heatmap_df <- sig_de %>%
  select(genes, logFC, adj.P.Val) %>%
  inner_join(term_gene_list, by = "genes")

# Optional: keep only strong changes
heatmap_df <- heatmap_df %>% filter(abs(logFC) > 1)

# Reshape to wide format: Terms as rows, genes as columns
heatmap_mat <- heatmap_df %>%
  select(Term, genes, logFC) %>%
  pivot_wider(names_from = genes, values_from = logFC, values_fill = 0)

# Convert back to long tidy format for ggplot heatmap
plot_df <- heatmap_mat %>%
  pivot_longer(-Term, names_to = "Gene", values_to = "logFC")

# Plot heatmap
p_heat <- ggplot(plot_df, aes(x = Gene, y = Term, fill = logFC)) +
  geom_tile(color = "grey70") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  labs(title = "Differential Expression across GO Terms",
       x = "Gene", y = "GO Term") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("Figures/GO_firing_terms_heatmap2.png", p_heat, width = 16, height = 4, dpi = 300)
