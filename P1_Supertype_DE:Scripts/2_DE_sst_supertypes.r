# run on conda activate de_env
setwd("P1_Supertype_DE")
library(ggrepel)
library(cowplot)
library(limma)
library(dplyr)
library(edgeR)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)


data <- readRDS("Files/Pseudobulk_SEAD_con.rds")

pseudobulk_metadata <- read.csv("Files/pseudobulk_metadata_SEAD_con.csv")%>%
  group_by(donor_id) %>%
  mutate(total_cells = sum(num_cells)) %>%
  ungroup() 

pseudobulk <- data %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID_Celltype") %>%
  separate(ID_Celltype, into = c("supertype", "donor_id"), sep = "_", extra = "merge", fill = "right") %>%
  filter(donor_id %in% pseudobulk_metadata$donor_id) %>%
  relocate(supertype, donor_id)

pseudobulk <- pseudobulk %>%
  mutate(supertype = gsub("-", "_", supertype))


CELL_TYPES <- c("Sst_1","Sst_4","Sst_5","Sst_7","Sst_3",
                "Sst_19","Sst_9","Sst_10","Sst_13","Sst_12",
                "Sst_11","Sst_20","Sst_22","Sst_23","Sst_25","Sst_2")

affected   <- c("Sst_3","Sst_20","Sst_22","Sst_25","Sst_2", "Sst_23", "Sst_11")
unaffected <- setdiff(CELL_TYPES, affected)

#collapse affeced vs unaffected sst supertypes
collapsed_counts <- pseudobulk %>%
  filter(supertype %in% CELL_TYPES) %>%
  mutate(Group = ifelse(supertype %in% affected, "Affected", "Unaffected")) %>%
  group_by(donor_id, Group) %>%
  summarise(across(-c(supertype), sum), .groups = "drop")

# Expression matrix
mat <- collapsed_counts %>%
  select(-donor_id, -Group) %>%
  t()

colnames(mat) <- paste(collapsed_counts$donor_id, collapsed_counts$Group, sep = "_")

# Metadata aligned
meta_collapsed <- collapsed_counts %>%
  select(donor_id, Group)
meta_collapsed$Group <- factor(meta_collapsed$Group, levels = c("Unaffected", "Affected"))

dge <- DGEList(counts = mat, genes = rownames(mat))

# Keep genes with ≥1 count in ≥80% of samples
min_samples <- ncol(mat) * 0.8
dge <- dge[rowSums(dge$counts >= 1) >= min_samples, ]

# Normalization
dge <- calcNormFactors(dge, method = "TMM")

#Model 
design <- model.matrix(~ Group, data = meta_collapsed)

vm <- voom(dge, design, plot = FALSE)
fit <- lmFit(vm, design)
fit <- eBayes(fit)

#Results

DE_aff_unaff <- topTable(
  fit,
  coef = "GroupAffected", 
  n = Inf,
  adjust.method = "BH"
)

# Save results
saveRDS(DE_aff_unaff, "Files/DE_results_Affected_vs_Unaffected.rds")

volcano <- EnhancedVolcano(
  DE_aff_unaff,
  lab = rownames(DE_aff_unaff),
  x = "logFC",
  y = "adj.P.Val",
  title = "Affected vs Unaffected Sst",
  pCutoff = 0.05,
  FCcutoff = 1
)

# Save to file
ggsave("Figures/Volcano_Affected_vs_Unaffected.png", plot = volcano,
       width = 8, height = 6, dpi = 300)


