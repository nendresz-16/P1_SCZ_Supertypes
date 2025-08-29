# conda activate crumblr_env
# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(ggplot2)
library(scattermore)
library(ggtree)
library(crumblr)
library(aplot)
library(tidyverse)
library(dreamlet)
library(kableExtra)
library(ggcorrplot)
library(RColorBrewer)
library(DelayedArray)
library(Seurat)
library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)
library(ggplot2)
library(ape)
setwd("P1_Compositional_analysis")



# ==== 1. Load & Filter Cohorts ====
filter_seurat_by_cells <- function(seu, min_cells = 500) {
  donor_counts <- table(seu$Donor)
  donors_to_keep <- names(donor_counts)[donor_counts >= min_cells]
  subset(seu, subset = Donor %in% donors_to_keep)
}

OFC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") %>%
  filter_seurat_by_cells()
Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") %>%
  filter_seurat_by_cells()
Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>%
  filter_seurat_by_cells()

# Split Ruz into McLean and MtSinai
Ruz_McLean  <- subset(Ruz, subset = Cohort == "McLean")
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")

# ==== 2. Keep Only Common Predicted IDs ====
ids_ofc    <- unique(OFC$predicted.id)
ids_bat    <- unique(Bat$predicted.id)
ids_mclean <- unique(Ruz_McLean$predicted.id)
ids_mtsinai<- unique(Ruz_MtSinai$predicted.id)

common_ids <- Reduce(intersect, list(ids_ofc, ids_bat, ids_mclean, ids_mtsinai))
cat("Shared predicted.ids across cohorts:\n")
print(common_ids)

OFC        <- subset(OFC, subset = predicted.id %in% common_ids)
Bat        <- subset(Bat, subset = predicted.id %in% common_ids)
Ruz_McLean <- subset(Ruz_McLean, subset = predicted.id %in% common_ids)
Ruz_MtSinai<- subset(Ruz_MtSinai, subset = predicted.id %in% common_ids)

# ==== 1. Subset objects by Age (<60 vs ≥60) ====
OFC_60         <- subset(OFC, subset = Age < 60)
Ruz_McLean_60  <- subset(Ruz_McLean, subset = Age < 60)
Ruz_MtSinai_60 <- subset(Ruz_MtSinai, subset = Age < 60)

OFC_60p        <- subset(OFC, subset = Age >= 60)
McLean_60p     <- subset(Ruz_McLean, subset = Age >= 60)
MtSinai_60p    <- subset(Ruz_MtSinai, subset = Age >= 60)

# Create lists for analysis
cohort_list_60 <- list(
  OFC        = OFC_60,
  McLean     = Ruz_McLean_60,
  MtSinai    = Ruz_MtSinai_60
)

cohort_list_60p <- list(
  OFC        = OFC_60p,
  McLean     = McLean_60p,
  MtSinai    = MtSinai_60p
)

# ==== 2. Define the dream formula ====
dream_formula <- ~ Diagnosis + Sex + scale(Age) + scale(PMI)

# ==== 3. Crumblr runner ====
run_crumblr_per_cohort <- function(seu, cohort_name, age_group) {
  message("→ Running Crumblr on ", cohort_name, " | Age group: ", age_group)

  # Convert Seurat -> SingleCellExperiment
  sce <- as.SingleCellExperiment(seu)

  # Aggregate to pseudo-bulk
  pb <- aggregateToPseudoBulk(
    sce,
    assay = "counts",
    cluster_id = "predicted.id", # adjust if renamed
    sample_id  = "Donor"
  )

  # Counts & metadata
  counts <- cellCounts(pb)
  meta   <- as.data.frame(colData(pb))

  # Ensure meta rows match counts
  meta <- meta[rownames(counts), , drop = FALSE]

  # Run Crumblr
  cobj <- crumblr(counts)
  fit <- dream(cobj, dream_formula, meta)
  fit <- eBayes(fit)

  # Extract results
  res <- topTable(fit, coef = "scale(Age)", number = Inf) %>%
    dplyr::select(logFC, AveExpr, t, P.Value, adj.P.Val) %>%
    tibble::rownames_to_column("CellType")

  # Add labels
  res$Cohort <- cohort_name
  res$AgeGroup <- age_group

  # Save per cohort + age group
  out_path <- paste0("Files/crumbler_results_", cohort_name, "_", age_group, ".csv")
  write_csv(res, out_path)

  message("✓ Saved: ", out_path)
  return(res)
}

# ==== 4. Run Crumblr for both <60 and ≥60 ====

# For <60 group
results_60 <- purrr::map2_dfr(
  .x = cohort_list_60,
  .y = names(cohort_list_60),
  .f = ~ run_crumblr_per_cohort(.x, .y, "Under60")
)

# For ≥60 group
results_60p <- purrr::map2_dfr(
  .x = cohort_list_60p,
  .y = names(cohort_list_60p),
  .f = ~ run_crumblr_per_cohort(.x, .y, "Over60")
)


# ============================
# Meta-analysis of Crumblr outputs by CellType
# ============================

# Activate conda env before starting
# conda activate r_env_meta_analysis

setwd("P1_Compositional_analysis")

library(dplyr)
library(metafor)
library(dplyr)
library(ggplot2)
library(patchwork)


# ============================
# 2. Read ONLY age-stratified files
# ============================
files <- list.files(
  path = "Files",
  pattern = "^crumbler_results_.*_(Under60|Over60)\\.csv$",
  full.names = TRUE
)

# Read and parse cohort + age group
all_data <- do.call(rbind, lapply(files, function(f) {
  df <- read.csv(f)
  
  # Extract cohort name from filename
  df$cohort <- sub("^crumbler_results_(.*?)_(Under60|Over60)\\.csv$", "\\1", basename(f))
  
  # Extract age group
  df$AgeGroup <- sub("^crumbler_results_.*?_(Under60|Over60)\\.csv$", "\\1", basename(f))
  
  # Calculate SE
  df$SE <- df$logFC / df$t
  df
}))

# ============================
# 3. Check structure
# ============================
cat("Total rows:", nrow(all_data), "\n")
cat("Cohorts included:", paste(unique(all_data$cohort), collapse = ", "), "\n")
cat("Age groups included:", paste(unique(all_data$AgeGroup), collapse = ", "), "\n")

# ============================
# 4. Run meta-analysis separately by age group
# ============================
meta_results <- list()

for (age_group in unique(all_data$AgeGroup)) {
  
  message("Running meta-analysis for Age Group: ", age_group)
  
  df_age <- subset(all_data, AgeGroup == age_group)
  cell_types <- unique(df_age$CellType)
  
  for (ct in cell_types) {
    df <- subset(df_age, CellType == ct & !is.na(SE) & SE != 0)
    
    # Require ≥2 cohorts per cell type
    if (nrow(df) >= 2) {
      res <- tryCatch({
        model <- rma(
          yi = df$logFC,
          sei = df$SE,
          method = "FE"
        )
        
        data.frame(
          CellType = ct,
          AgeGroup = age_group,
          estimate = model$b,
          se = model$se,
          zval = model$zval,
          pval = model$pval,
          ci.lb = model$ci.lb,
          ci.ub = model$ci.ub,
          k = model$k,
          I2 = model$I2
        )
      }, error = function(e) NULL)
      
      if (!is.null(res)) {
        meta_results[[paste(ct, age_group, sep = "_")]] <- res
      }
    }
  }
}

# ============================
# 5. Combine results
# ============================
final_results <- do.call(rbind, meta_results)
# Split final results into Under60 and Over60
final_results_under60 <- subset(final_results, AgeGroup == "Under60")
final_results_over60  <- subset(final_results, AgeGroup == "Over60")

# Adjust p-values separately for each group
final_results_under60$padj <- p.adjust(final_results_under60$pval, method = "fdr")
final_results_over60$padj  <- p.adjust(final_results_over60$pval, method = "fdr")


# Under60 significant results
sig_under60 <- subset(final_results_under60, padj < 0.1)
cat("\nSignificant cell types (Under60, FDR < 0.1):\n")
print(unique(sig_under60$CellType))

# Over60 significant results
sig_over60 <- subset(final_results_over60, padj < 0.1)
cat("\nSignificant cell types (Over60, FDR < 0.1):\n")
print(unique(sig_over60$CellType))


# Keep only Sst subtypes for each group
sst_under60 <- final_results_under60 %>%
  filter(grepl("^Sst", CellType))

sst_over60 <- final_results_over60 %>%
  filter(grepl("^Sst", CellType))

sst_combined <- bind_rows(
  sst_under60 %>% mutate(AgeGroup = "Under60"),
  sst_over60 %>% mutate(AgeGroup = "Over60")
)



# Define custom Sst order
sst_order <- c(
  "Sst_1", "Sst_4", "Sst_5", "Sst_7",
  "Sst_3", "Sst_19", "Sst_9", "Sst_10",
  "Sst_13", "Sst_12", "Sst_11", "Sst_20",
  "Sst_22", "Sst_23", "Sst_25", "Sst_2"
)

# ----------------------------
# 1. Subset Under60 results
# ----------------------------
sst_under60 <- final_results_under60 %>%
  filter(CellType %in% sst_order) %>%
  mutate(
    CellType = factor(CellType, levels = rev(sst_order)),
    sig = case_when(
      padj < 0.01 ~ "***",
      padj < 0.05 ~ "**",
      padj < 0.1 ~ "*",
      TRUE ~ ""
    )
  )

# ----------------------------
# 2. Subset Over60 results
# ----------------------------
sst_over60 <- final_results_over60 %>%
  filter(CellType %in% sst_order) %>%
  mutate(
    CellType = factor(CellType, levels = rev(sst_order)),
    sig = case_when(
      padj < 0.01 ~ "***",
      padj < 0.05 ~ "**",
      padj < 0.1 ~ "*",
      TRUE ~ ""
    )
  )

# ----------------------------
# 3. Plot for Under60
# ----------------------------
p1 <- ggplot(sst_under60, aes(x = CellType, y = estimate)) +
  geom_col(fill = "#a89cf6ff") +
  geom_text(aes(label = sig, vjust = ifelse(estimate >= 0, -0.5, 1.5)), size = 5) +
  labs(
    x = "Sst Subtype",
    y = "Meta-analysis Estimate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------
# 4. Plot for Over60
# ----------------------------
p2 <- ggplot(sst_over60, aes(x = CellType, y = estimate)) +
  geom_col(fill = "#ecf085ff") +
  geom_text(aes(label = sig, vjust = ifelse(estimate >= 0, -0.5, 1.5)), size = 5) +
  labs(
    x = "Sst Subtype",
    y = "Meta-analysis Estimate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------
# 5. Print both plots
# ----------------------------

library(patchwork)  

# Combine plots vertically
combined_plot <- p1 / p2   

# Save combined figure
ggsave("Figures/Sst_meta_analysis_age_stratified.png", combined_plot, width = 12, height = 10)