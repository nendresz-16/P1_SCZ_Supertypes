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




# ==== 1. Load SEAD and Donor ≥ 500 cells ====

SEAD <- readRDS("/project/rrg-shreejoy/nendresz/SEAD_data/SEAD_obj_gene_symbols.rds") 

# Rename columns
SEAD@meta.data <- SEAD@meta.data %>%
  rename(
    Donor = donor_id,
    Age   = `Age.at.death`,
    Sex   = sex
  )

# Recode Cognitive.status values
SEAD@meta.data$Diagnosis <- recode(SEAD@meta.data$Cognitive.status,
                                   "Dementia" = "AD",
                                   "No dementia" = "Control")

donor_counts <- table(SEAD$Donor)
donors_to_keep <- names(donor_counts)[donor_counts >= 500]
SEAD <- subset(SEAD, subset = Donor %in% donors_to_keep)

# Combine all counts layers into one
mat_list <- lapply(Layers(SEAD[["RNA"]]), function(layer) {
  LayerData(SEAD[["RNA"]], layer = layer)
})
mat <- do.call(cbind, mat_list)

# Fix rownames globally
rownames(mat) <- make.unique(rownames(mat))

# Recreate a clean assay
SEAD[["RNA"]] <- CreateAssayObject(counts = mat)
DefaultAssay(SEAD) <- "RNA"
Layers(SEAD[["RNA"]])



# ==== 2. Define model formula ====
dream_formula <- ~ Diagnosis + Sex + scale(Age) + scale(PMI)

# ==== 3. Run Crumblr ====

sce <- as.SingleCellExperiment(SEAD)

pb <- aggregateToPseudoBulk(
  sce,
  assay = "counts",
  cluster_id = "Supertype",
  sample_id  = "Donor"
)

counts <- cellCounts(pb)
meta   <- as.data.frame(colData(pb))
meta   <- meta[rownames(counts), , drop = FALSE]

cobj <- crumblr(counts)
fit <- dream(cobj, dream_formula, meta)
fit <- eBayes(fit)

res <- topTable(fit, coef = "DiagnosisAD", number = Inf) %>%
  dplyr::select(logFC, AveExpr, t, P.Value, adj.P.Val) %>%
  tibble::rownames_to_column("CellType")

out_path <- "Files/crumbler_results_SEAD.csv"
write_csv(res, out_path)

