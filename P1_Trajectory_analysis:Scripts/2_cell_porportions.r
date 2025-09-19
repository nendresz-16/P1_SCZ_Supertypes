# conda activate r_anndata_env
# -------------------------------------------------------- #
# Packages
# -------------------------------------------------------- #
library(dplyr)
library(tibble)
library(reshape2)   
Sys.setenv(RETICULATE_PYTHON = "/scratch/nendresz/envs/phate_env/bin/python")
library(reticulate)
library(anndata) 
setwd("P1_Trajectory_analysis")
# -------------------------------------------------------- #
# Compute subpopulation proportions                        #
# -------------------------------------------------------- #
Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds")

# Cell counts per Donor x predicted.id + donor-wise prevalence
df <- Bat@meta.data %>%
  count(predicted.id, Donor) %>%
  group_by(Donor) %>%
  mutate(prevalence = n / sum(n)) %>%
  ungroup()
gc()

# Wide matrices: Donor x predicted.id
proportions <- dcast(
  df, Donor ~ predicted.id,
  value.var = "prevalence", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor") %>% as.matrix()

counts <- dcast(
  df, Donor ~ predicted.id,
  value.var = "n", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor") %>% as.matrix()

saveRDS(proportions, "Files/subpopulation.proportion.matrix.rds")
saveRDS(counts,      "Files/subpopulation.counts.matrix.rds")

# -------------------------------------------------------- #
# Create base AnnData (participants x subpopulations)
# -------------------------------------------------------- #
ids  <- rownames(proportions)          # donors
feat <- colnames(proportions)          # predicted.id labels

stopifnot(all(rownames(proportions) == rownames(counts)))
stopifnot(all(colnames(proportions) == colnames(counts)))

# Minimal var/obs
var_tbl <- data.frame(row.names = feat)
obs_tbl <- data.frame(row.names = ids)

data <- AnnData(
  X      = proportions,
  obs    = obs_tbl,
  var    = var_tbl,
  layers = list(
    counts    = counts[ids, feat, drop = FALSE],
    sqrt.prev = sqrt(proportions)
  )
)

anndata::write_h5ad(data, "Files/subpopulation.proportions.h5ad")

# -------------------------------------------------------- #
# Append participant metadata (donor-level, aligned)       #
# -------------------------------------------------------- #
# Build a donor-level table (e.g., # cells per donor) so indexing matches donors
donor_meta <- Bat@meta.data %>%
  group_by(Donor) %>%
  summarise(
    Sex = unique(Sex)[1],
    Age = unique(Age)[1],
    Diagnosis = unique(Diagnosis)[1],
    .groups = "drop"
  ) %>%
  left_join(
    Bat@meta.data %>%
      count(Donor, name = "n_cells"),
    by = "Donor"
  ) %>%
  column_to_rownames("Donor")

data <- anndata::read_h5ad("Files/subpopulation.proportions.h5ad")
donor_meta <- donor_meta[data$obs_names, , drop = FALSE]   # align to donors
data$obs[names(donor_meta)] <- donor_meta
data$obsm$meta.data <- donor_meta
anndata::write_h5ad(data, "Files/subpopulation.proportions.h5ad")




rm(donor_meta)

# -------------------------------------------------------- #
# Cell type level of aggregation (store in data$uns)       #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("Files/subpopulation.proportions.h5ad")

cell.type.df <- df %>%
  group_by(Donor, predicted.id) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Donor) %>%
  mutate(prevalence = n / sum(n)) %>%
  ungroup()

ct.counts <- dcast(
  cell.type.df, Donor ~ predicted.id,
  value.var = "n", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor")

ct.prev <- dcast(
  cell.type.df, Donor ~ predicted.id,
  value.var = "prevalence", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor")

# Align to AnnData (row order = donors, col order = X columns)
ct.counts <- ct.counts[data$obs_names, colnames(data$X), drop = FALSE]
ct.prev   <- ct.prev[  data$obs_names, colnames(data$X), drop = FALSE]

data$uns$cell.types <- list(
  counts    = as.matrix(ct.counts),
  prev      = as.matrix(ct.prev),
  sqrt.prev = sqrt(as.matrix(ct.prev))
)

anndata::write_h5ad(data, "Files/subpopulation.proportions.h5ad")


##########################################
##########################################
############RUZICKA####################
########################################




# conda activate r_anndata_env
# -------------------------------------------------------- #
# Packages
# -------------------------------------------------------- #
library(dplyr)
library(tibble)
library(reshape2)   
Sys.setenv(RETICULATE_PYTHON = "/scratch/nendresz/envs/phate_env/bin/python")
library(reticulate)
library(anndata) 
setwd("P1_Trajectory_analysis")
# -------------------------------------------------------- #
# Compute subpopulation proportions                        #
# -------------------------------------------------------- #
Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 

# Cell counts per Donor x predicted.id + donor-wise prevalence
df <- Ruz@meta.data %>%
  count(predicted.id, Donor) %>%
  group_by(Donor) %>%
  mutate(prevalence = n / sum(n)) %>%
  ungroup()
gc()

# Wide matrices: Donor x predicted.id
proportions <- dcast(
  df, Donor ~ predicted.id,
  value.var = "prevalence", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor")  %>%  mutate(across(everything(), as.numeric)) %>%  
  as.matrix()

counts <- dcast(
  df, Donor ~ predicted.id,
  value.var = "n", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor") %>% as.matrix()

saveRDS(proportions, "Files/subpopulation.proportion.matrix.rds")
saveRDS(counts,      "Files/subpopulation.counts.matrix.rds")

# -------------------------------------------------------- #
# Create base AnnData (participants x subpopulations)
# -------------------------------------------------------- #
ids  <- rownames(proportions)          # donors
feat <- colnames(proportions)          # predicted.id labels

stopifnot(all(rownames(proportions) == rownames(counts)))
stopifnot(all(colnames(proportions) == colnames(counts)))

# Minimal var/obs
var_tbl <- data.frame(row.names = feat)
obs_tbl <- data.frame(row.names = ids)

data <- AnnData(
  X      = proportions,
  obs    = obs_tbl,
  var    = var_tbl,
  layers = list(
    counts    = counts[ids, feat, drop = FALSE],
    sqrt.prev = sqrt(proportions)
  )
)

anndata::write_h5ad(data, "Files/subpopulation.proportions_Ruz.h5ad")

# -------------------------------------------------------- #
# Append participant metadata (donor-level, aligned)       #
# -------------------------------------------------------- #
# Build a donor-level table (e.g., # cells per donor) so indexing matches donors
donor_meta <- Ruz@meta.data %>%
  group_by(Donor) %>%
  summarise(
    Sex = unique(Sex)[1],
    Age = unique(Age)[1],
    cohort = unique(cohort)[1],
    Diagnosis = unique(Diagnosis)[1],
    .groups = "drop"
  ) %>%
  left_join(
    Ruz@meta.data %>%
      count(Donor, name = "n_cells"),
    by = "Donor"
  ) %>%
  column_to_rownames("Donor")

data <- anndata::read_h5ad("Files/subpopulation.proportions_Ruz.h5ad")
donor_meta <- donor_meta[data$obs_names, , drop = FALSE]   # align to donors
data$obs[names(donor_meta)] <- donor_meta
data$obsm$meta.data <- donor_meta
anndata::write_h5ad(data, "Files/subpopulation.proportions_Ruz.h5ad")




rm(donor_meta)

# -------------------------------------------------------- #
# Cell type level of aggregation (store in data$uns)       #
# -------------------------------------------------------- #
data <- anndata::read_h5ad("Files/subpopulation.proportions_Ruz.h5ad")

cell.type.df <- df %>%
  group_by(Donor, predicted.id) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Donor) %>%
  mutate(prevalence = n / sum(n)) %>%
  ungroup()

ct.counts <- dcast(
  cell.type.df, Donor ~ predicted.id,
  value.var = "n", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor")

ct.prev <- dcast(
  cell.type.df, Donor ~ predicted.id,
  value.var = "prevalence", fill = 0, fun.aggregate = sum
) %>% tibble::column_to_rownames("Donor")

# Align to AnnData (row order = donors, col order = X columns)
ct.counts <- ct.counts[data$obs_names, colnames(data$X), drop = FALSE]
ct.prev   <- ct.prev[  data$obs_names, colnames(data$X), drop = FALSE]

data$uns$cell.types <- list(
  counts    = as.matrix(ct.counts),
  prev      = as.matrix(ct.prev),
  sqrt.prev = sqrt(as.matrix(ct.prev))
)

anndata::write_h5ad(data, "Files/subpopulation.proportions_Ruz.h5ad")
