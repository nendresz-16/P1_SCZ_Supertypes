library(Seurat) 
library(dplyr)
library(ggplot2)
library(tibble)
setwd("P1_Supertype_DE")

counts <- readRDS("/project/rrg-shreejoy/nendresz/raw_counts_ref.rds")

meta <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")

obj <- CreateSeuratObject(counts = t(counts), meta.data = meta)

levels(factor(obj@meta.data$Cognitive.status))




Idents(obj) <- "Supertype"


CELL_TYPES <- c("Sst_1","Sst_4","Sst_5","Sst_7","Sst_3",
                "Sst_19","Sst_9","Sst_10","Sst_13","Sst_12",
                "Sst_11","Sst_20","Sst_22","Sst_23","Sst_25","Sst_2")



obj_sst <- subset(obj, idents = CELL_TYPES)

obj_sst <- NormalizeData(obj_sst, assay = "RNA")
obj_sst <- FindVariableFeatures(obj_sst, assay = "RNA")
obj_sst <- ScaleData(obj_sst, assay = "RNA")

markers_sst <- FindAllMarkers(
  obj_sst,
  only.pos = TRUE,          # only return positive markers
  min.pct = 0.25,           # expressed in at least 25% of cells
  logfc.threshold = 0.25    # minimum log-fold change
)

# Inspect the results
head(markers_sst)

top10 <- markers_sst %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()


# Your preferred order
sst_order <- c(
  "Sst_1", "Sst_4", "Sst_5", "Sst_7",
  "Sst_3", "Sst_19", "Sst_9", "Sst_10",
  "Sst_13", "Sst_12", "Sst_11", "Sst_20",
  "Sst_22", "Sst_23", "Sst_25", "Sst_2"
)

# Make sure identities follow this order
Idents(obj_sst) <- factor(Idents(obj_sst), levels = sst_order)

# Now DoHeatmap will respect this order
p <- DoHeatmap(obj_sst, features = top10$gene, size = 3) + NoLegend()


# Save as PNG
ggsave("Figures/Sst_marker_heatmap.png", p, width = 16, height = 18, dpi = 300)







#########################################################
counts <- readRDS("/project/rrg-shreejoy/nendresz/raw_counts_ref.rds")

meta <- readRDS("/project/rrg-shreejoy/nendresz/Neurotypical_ref_metadata.rds")

obj <- CreateSeuratObject(counts = t(counts), meta.data = meta)


Idents(obj) <- "Supertype"


CELL_TYPES <- c("Sst_1","Sst_4","Sst_5","Sst_7","Sst_3",
                "Sst_19","Sst_9","Sst_10","Sst_13","Sst_12",
                "Sst_11","Sst_20","Sst_22","Sst_23","Sst_25","Sst_2")



obj_sst <- subset(obj, idents = CELL_TYPES)

obj_sst <- NormalizeData(obj_sst, assay = "RNA")
obj_sst <- FindVariableFeatures(obj_sst, assay = "RNA")
obj_sst <- ScaleData(obj_sst, assay = "RNA")

# Affected vs unaffected

# Define affected supertypes
affected_ssts <- c("Sst_2","Sst_3","Sst_22","Sst_23","Sst_25","Sst_11")

# Add a new column: affected vs unaffected
obj_sst$Affected_group <- ifelse(
  obj_sst$Supertype %in% affected_ssts, 
  "Affected", 
  "Unaffected"
)

# Set identities to this new grouping
Idents(obj_sst) <- "Affected_group"

# Check group sizes
table(Idents(obj_sst))


markers_aff_vs_unaff <- FindMarkers(
  obj_sst,
  ident.1 = "Affected",
  ident.2 = "Unaffected",
  assay = "RNA",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(markers_aff_vs_unaff)

# Take top 10 markers by logFC
top10 <- markers_aff_vs_unaff %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10)
top10 <- top10 %>% rownames_to_column("gene")

# Plot
p1 <- DoHeatmap(obj_sst, features = top10$gene, size = 3) + theme(legend.position = "right")

ggsave("Figures/Affected_vs_Unaffected_Sst_markers.png", p1, width = 12, height = 8, dpi = 300) 
