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




HBCC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/HBCC_updated.rds")
MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
# Merge

MSSM <- merge(MSSM1, y = MSSM2)
MSSM <- JoinLayers(MSSM)

Layers(MSSM[["RNA"]])



# ==== 1. Load and subset cohorts ====
OFC   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") 
Bat   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") 
Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 


# Split Ruz into McLean and MtSinai
Ruz_McLean  <- subset(Ruz, subset = Cohort == "McLean")
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


# ==== 1. Load and subset cohorts (Age < 70 and Donor ≥ 500 cells) ====

filter_seurat_by_age_and_cells <- function(seu) {
  seu <- subset(
    seu, 
    subset = Age < 70 & grepl("^Sst_", predicted.id)
  )
  return(seu)
}

OFC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") %>%
  filter_seurat_by_age_and_cells()

Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") %>%
  filter_seurat_by_age_and_cells()

Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>%
  filter_seurat_by_age_and_cells()

# Split Ruz into McLean and MtSinai
Ruz_McLean  <- subset(Ruz, subset = Cohort == "McLean")
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


# === 1. Get unique predicted.ids in each object ===
ids_ofc   <- unique(OFC$predicted.id)
ids_bat   <- unique(Bat$predicted.id)
ids_mclean <- unique(Ruz_McLean$predicted.id)
ids_mtsinai <- unique(Ruz_MtSinai$predicted.id)
ids_hbcc <- unique(HBCC$predicted.id)
ids_mssm <- unique(MSSM$predicted.id)

# === 2. Find common predicted.ids across all cohorts ===
common_ids <- Reduce(intersect, list(ids_ofc, ids_bat, ids_mclean, ids_mtsinai, ids_hbcc, ids_mssm))

# Optional: print them
cat("Shared predicted.ids across cohorts:\n")
print(common_ids)

# === 3. Subset each Seurat object ===
OFC <- subset(OFC, subset = predicted.id %in% common_ids)
Bat <- subset(Bat, subset = predicted.id %in% common_ids)
Ruz_McLean <- subset(Ruz_McLean, subset = predicted.id %in% common_ids)
Ruz_MtSinai <- subset(Ruz_MtSinai, subset = predicted.id %in% common_ids)
HBCC <- subset(HBCC, subset = predicted.id %in% common_ids)
MSSM <- subset(MSSM, subset = predicted.id %in% common_ids)


# Put in a named list
cohort_list <- list(
  OFC        = OFC,
  Bat        = Bat,
  McLean     = Ruz_McLean,
  MtSinai    = Ruz_MtSinai,
  HBCC       = HBCC,
  MSSSM      = MSSM
)

# ==== 2. Define the model formula ====
dream_formula <- ~ Diagnosis + Sex + scale(Age) + scale(PMI)

# ==== 3. Function to run Crumblr and export results per cohort ====
run_crumblr_per_cohort <- function(seu, cohort_name) {
  message("→ Running Crumblr on ", cohort_name)

  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seu)

  # Aggregate to pseudo-bulk
  pb <- aggregateToPseudoBulk(
    sce,
    assay = "counts",
    cluster_id = "predicted.id",  # change to "Subtype" if renamed
    sample_id  = "Donor"
  )

  # Get counts and metadata
  counts <- cellCounts(pb)
  meta   <- as.data.frame(colData(pb))

  # Run Crumblr and dream model
  cobj <- crumblr(counts)
  meta <- meta[rownames(counts), , drop = FALSE]

  # Fit model
  fit <- dream(cobj, dream_formula, meta)
  fit <- eBayes(fit)

  # Get results and export
  res <- topTable(fit, coef = "DiagnosisSchizophrenia", number = Inf) %>%
    dplyr::select(logFC, AveExpr, t, P.Value, adj.P.Val)%>%
    tibble::rownames_to_column("CellType")

  # Save to CSV
  out_path <- paste0("Files/crumbler_results_SST", cohort_name, ".csv")
  write_csv(res, out_path)
  message("✓ Saved results to ", out_path)

  return(invisible(res))
}

# ==== 4. Run across all cohorts ====
purrr::walk2(
  .x = cohort_list,
  .y = names(cohort_list),
  .f = run_crumblr_per_cohort
)



# Function to extract donor count
get_sample_size <- function(seurat_obj, donor_col = "Donor") {
  if (donor_col %in% colnames(seurat_obj@meta.data)) {
    return(length(unique(seurat_obj@meta.data[[donor_col]])))
  } else {
    warning(paste("Column", donor_col, "not found in metadata"))
    return(NA)
  }
}

# Get sample sizes
sample_sizes <- list(
  OFC        = get_sample_size(OFC),
  Bat        = get_sample_size(Bat),
  Ruz_McLean = get_sample_size(Ruz_McLean),
  Ruz_MtSinai = get_sample_size(Ruz_MtSinai),
  HBCC        = get_sample_size(HBCC),
  MSSM        = get_sample_size(MSSM)
)

# Print nicely
print(sample_sizes)



# conda activate r_env_meta_analysis
setwd("P1_Compositional_analysis")
library(dplyr)
# List your CSV files
files <- c(
  "Files/crumbler_results_SSTBat.csv",
  "Files/crumbler_results_SSTHBCC.csv",
  "Files/crumbler_results_SSTMcLean.csv",
  "Files/crumbler_results_SSTMSSSM.csv",
  "Files/crumbler_results_SSTMtSinai.csv",
  "Files/crumbler_results_SSTOFC.csv"
)


# Read and combine
all_data <- do.call(rbind, lapply(files, function(f) {
  df <- read.csv(f)
  df$cohort <- sub("crumbler_results_(.*)\\.csv", "\\1", f)
  df$SE <- df$logFC / df$t  # compute standard error
  return(df)
}))


library(metafor)

# Get unique cell types
cell_types <- unique(all_data$CellType)

# Empty list to store results
meta_results <- list()

# Loop through each cell type
for (ct in cell_types) {
  df <- subset(all_data, CellType == ct & !is.na(SE) & SE != 0)

  if (nrow(df) >= 2) {  # need at least 2 cohorts for meta-analysis
    res <- tryCatch({
      model <- rma(yi = df$logFC, sei = df$SE, method = "FE")
      data.frame(CellType = ct,
                 estimate = model$b,
                 se = model$se,
                 zval = model$zval,
                 pval = model$pval,
                 ci.lb = model$ci.lb,
                 ci.ub = model$ci.ub)
    }, error = function(e) NULL)

    if (!is.null(res)) {
      meta_results[[ct]] <- res
    }
  }
}

# Combine results into a single data.frame
final_results <- do.call(rbind, meta_results)



final_results$padj <- p.adjust(final_results$pval, method = "fdr")



sig_results <- subset(final_results, padj < 0.1)


unique(sig_results$CellType)

# === 3. Prepare meta-analysis for plotting ===
final_results$Cohort <- "Meta-analysis"
final_results$padj <- p.adjust(final_results$pval, method = "fdr")
final_results$neglogP <- -log10(final_results$pval)
final_results$sig_label <- dplyr::case_when(
  final_results$padj < 0.05 ~ "+",
  final_results$padj < 0.1  ~ "/",
  TRUE ~ ""
)
meta_plot_data <- final_results[, c("CellType", "estimate", "padj", "neglogP", "sig_label", "Cohort")]

# === 4. Prepare cohort-level data ===
all_data$padj <- as.numeric(all_data$adj.P.Val)
all_data$neglogP <- -log10(all_data$padj)
all_data$sig_label <- dplyr::case_when(
  all_data$padj < 0.05 ~ "+",
  all_data$padj < 0.1  ~ "/",
  TRUE ~ ""
)
colnames(all_data)[colnames(all_data) == "cohort"] <- "Cohort"

cohort_plot_data <- all_data[, c("CellType", "logFC", "padj", "neglogP", "sig_label", "Cohort")]
colnames(cohort_plot_data)[colnames(cohort_plot_data) == "logFC"] <- "estimate"

# === 5. Combine meta + cohorts ===
plot_data <- rbind(meta_plot_data, cohort_plot_data)
plot_data$Cohort <- factor(plot_data$Cohort, levels = c("Meta-analysis", sort(unique(cohort_plot_data$Cohort))))
plot_data$CellType <- factor(plot_data$CellType, levels = rev(sort(unique(plot_data$CellType))))
plot_data$Cohort <- recode(plot_data$Cohort, "Bat" = "Batiuk","OFC" = "Fröhlich")



library(ggplot2)

# === 6. Plot ===
p <- ggplot(plot_data, aes(x = CellType, y = Cohort, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.3) +  # white grid lines
  geom_text(aes(label = sig_label), size = 4) +  # significance "+" on top
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "Effect Size"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Cell Type",
    y = "Cohort"
  )

# === 7. Save plot ===
ggsave("Figures/SCZ_meta_cohort_heatmap_in_SST.png", plot = p, width = 8, height = 3, dpi = 300)






