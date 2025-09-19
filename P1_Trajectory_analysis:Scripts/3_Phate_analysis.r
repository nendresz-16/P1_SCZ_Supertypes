#conda activate /scratch/nendresz/envs/phate_env
#module load r/4.3.1  
setwd("P1_Trajectory_analysis")

# Point R to your conda Python
Sys.setenv(RETICULATE_PYTHON = "/scratch/nendresz/envs/phate_env/bin/python")

library(reticulate)
library(dplyr)

py_config()

sc <- import("scanpy")
ad <- import("anndata")
ph <- import("phate")

data <- ad$read_h5ad("Files/subpopulation.proportions.h5ad")

# Participant clustering based on cellular environment representation
sc$pp$neighbors(data, n_neighbors = as.integer(2), use_rep = "X", metric = "cosine")
sc$tl$leiden(data, resolution =.25)
sc$tl$leiden(data, resolution = .75, restrict_to=reticulate::tuple("leiden", reticulate::np_array(c("0"))))
data$obs["clusters"] = plyr::mapvalues(data$obs$leiden_R, levels(data$obs$leiden_R), 1:length(levels(data$obs$leiden_R)))
data$obs$leiden <- data$obs$leiden_R <- NULL


data$obs$core <- !data$obs$clusters %in% c(9,10)

# 2D visualization of landscape
perp <- as.integer(min(10, floor((data$n_obs - 1)/3)))  # -> 7 when n_obs=23

# t-SNE (set perplexity; also lower learning rate for tiny n)
sc$tl$tsne(
  data,
  use_rep       = "X",
  perplexity    = perp,
  learning_rate = 50L
)


#sc$tl$tsne(data, n_pcs = 0, use_rep = "X", learning_rate = 100)
sc$tl$umap(data, maxiter = as.integer(3000), spread = 3)

# Save t-SNE plot colored by diagnosis
sc$pl$tsne(data, color="Diagnosis", show=FALSE)

# Save manually to your desired path
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Bat_tsne_diagnosis.png", dpi=300)

# Save UMAP plot colored by diagnosis
sc$pl$umap(data, color="Diagnosis", save="Figures/Bat_umap_diagnosis.png")


# Turn off Scanpy's autosave
sc$pl$umap(data, color="Diagnosis", show=FALSE)

# Save manually to your desired location
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Bat_umap_diagnosis.png", dpi=300)



sc$external$tl$phate(data, 
                     n_components = as.integer(3),  
                     k = as.integer(10), a = as.integer(40), 
                     knn_dist =  "euclidean", mds_dist = "euclidean", 
                     mds_solver = "smacof", verbose = F)

sc$pl$embedding(
  data,
  basis = "phate",           # Must match 'X_phate' in obsm
  color = "Diagnosis",
  projection = "3d",
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Bat_PHATE3D_Diagnosis.png", dpi = 300)


# 1. Run PHATE in 2D
sc$external$tl$phate(
  data,
  n_components = as.integer(2),
  knn = as.integer(10),
  decay = as.integer(40),
  knn_dist = "euclidean",
  mds_dist = "euclidean",
  mds_solver = "smacof",
  verbose = FALSE
)

# 2. Plot the 2D PHATE (this looks for data$obsm[["X_phate"]])
sc$pl$embedding(
  data,
  basis = "phate",         # This corresponds to X_phate
  color = "Diagnosis",     # Or "clusters", "core", etc.
  show = FALSE             # Don't display in R; we'll save manually
)

# 3. Save the plot manually using matplotlib
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Bat_PHATE2D_Diagnosis.png", dpi = 300)



sc$pl$embedding(
  data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Sst_22",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Bat_PHATE2D_Sst25.png", dpi = 300)



##################################################
##################################################
############                           ###########
############         RUZICKA           ###########
############                           ###########
##################################################
##################################################


#conda activate /scratch/nendresz/envs/phate_env
#module load r/4.3.1  
setwd("P1_Trajectory_analysis")

# Point R to your conda Python
Sys.setenv(RETICULATE_PYTHON = "/scratch/nendresz/envs/phate_env/bin/python")

library(reticulate)
library(dplyr)

py_config()

sc <- import("scanpy")
ad <- import("anndata")
ph <- import("phate")

data <- ad$read_h5ad("Files/subpopulation.proportions_Ruz.h5ad")

#Split into Mclean and Mtsinai 
obs_df <- reticulate::py_to_r(data$obs)
donor_names <- reticulate::py_to_r(data$obs_names)

# Get donor names (row names) by cohort
mt_sinai_donors <- donor_names[obs_df$cohort == "MtSinai"]
mclean_donors   <- donor_names[obs_df$cohort == "McLean"]

# Subset using donor IDs (not numeric indices)
mt_sinai_data <- data[mt_sinai_donors, ]
mclean_data   <- data[mclean_donors, ]

mt_sinai_data$write_h5ad("Files/subpopulation_Ruz_MtSinai.h5ad")
mclean_data$write_h5ad("Files/subpopulation_Ruz_McLean.h5ad")

library(reticulate)
sc  <- import("scanpy")
plt <- import("matplotlib.pyplot")
py_run_string("import numpy as np")  # needed for np_array
np <- import("numpy")


#------------------------------------------------------------#
# Function to process an AnnData object and save tSNE/UMAP
#------------------------------------------------------------#

# Load dependencies
ad <- import("anndata")
plt <- import("matplotlib.pyplot")
np  <- import("numpy")
sc  <- import("scanpy")

# Define plotting + clustering function
process_and_plot <- function(data, tag = "default") {
  sc$pp$neighbors(data, n_neighbors = as.integer(2), use_rep = "X", metric = "cosine")
  sc$tl$leiden(data, resolution = 0.25)
  sc$tl$leiden(
    data,
    resolution = 0.75,
    restrict_to = reticulate::tuple("leiden", np$array(c("0")))
  )

  data$obs["clusters"] <- plyr::mapvalues(
    data$obs$leiden_R,
    levels(data$obs$leiden_R),
    1:length(levels(data$obs$leiden_R))
  )
  data$obs$leiden <- NULL
  data$obs$leiden_R <- NULL

  data$obs$core <- !data$obs$clusters %in% c(9, 10)

  perp <- as.integer(min(10, floor((data$n_obs - 1)/3)))
  sc$tl$tsne(data, use_rep = "X", perplexity = perp, learning_rate = 50L)
  sc$tl$umap(data, maxiter = as.integer(3000), spread = 3)

  sc$pl$tsne(data, color = "Diagnosis", show = FALSE)
  plt$savefig(paste0("Figures/", tag, "_tsne_diagnosis.png"), dpi = 300)

  sc$pl$umap(data, color = "Diagnosis", show = FALSE)
  plt$savefig(paste0("Figures/", tag, "_umap_diagnosis.png"), dpi = 300)
}

# ------------------------------
# Load and assign to variables
# ------------------------------

mclean_data   <- ad$read_h5ad("Files/subpopulation_Ruz_McLean.h5ad")
mt_sinai_data <- ad$read_h5ad("Files/subpopulation_Ruz_MtSinai.h5ad")

# Optional: Run PHATE or clustering/plotting for each
process_and_plot(mclean_data,   tag = "Ruz_McLean")
process_and_plot(mt_sinai_data, tag = "Ruz_MtSinai")



sc$external$tl$phate(mclean_data, 
                     n_components = as.integer(3),  
                     k = as.integer(10), a = as.integer(40), 
                     knn_dist =  "euclidean", mds_dist = "euclidean", 
                     mds_solver = "smacof", verbose = F)

sc$pl$embedding(mclean_data,  basis = "phate", color = "Diagnosis", projection = "3d", show = FALSE)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mclean_PHATE3D_Diagnosis.png", dpi = 300)


# 1. Run PHATE in 2D
sc$external$tl$phate(mclean_data, n_components = as.integer(2), knn = as.integer(10), decay = as.integer(40), 
  knn_dist = "euclidean",
  mds_dist = "euclidean",
  mds_solver = "smacof",
  verbose = FALSE
)

# 2. Plot the 2D PHATE (this looks for data$obsm[["X_phate"]])
sc$pl$embedding(
  mclean_data,
  basis = "phate",         # This corresponds to X_phate
  color = "Diagnosis",     # Or "clusters", "core", etc.
  show = FALSE             # Don't display in R; we'll save manually
)

# 3. Save the plot manually using matplotlib
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mclean_PHATE2D_Diagnosis.png", dpi = 300)



sc$pl$embedding(
   mclean_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Sst_25",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mclean_PHATE2D_Sst25.png", dpi = 300)


sc$pl$embedding(
    mclean_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Age",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mclean_PHATE2D_Age.png", dpi = 300)

sc$pl$embedding(
    mclean_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Astro_2",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mclean_PHATE2D_Astro2.png", dpi = 300)

sc$pl$embedding(
    mclean_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Micro-PVM_2",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mclean_PHATE2D_Micro2
.png", dpi = 300)






#########################################################



sc$external$tl$phate(mt_sinai_data, 
                     n_components = as.integer(3),  
                     k = as.integer(10), a = as.integer(40), 
                     knn_dist =  "euclidean", mds_dist = "euclidean", 
                     mds_solver = "smacof", verbose = F)

sc$pl$embedding(mt_sinai_data,  basis = "phate", color = "Diagnosis", projection = "3d", show = FALSE)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE3D_Diagnosis.png", dpi = 300)


# 1. Run PHATE in 2D
sc$external$tl$phate(mt_sinai_data, n_components = as.integer(2), knn = as.integer(10), decay = as.integer(40), 
  knn_dist = "euclidean",
  mds_dist = "euclidean",
  mds_solver = "smacof",
  verbose = FALSE
)

# 2. Plot the 2D PHATE (this looks for data$obsm[["X_phate"]])
sc$pl$embedding(
  mt_sinai_data,
  basis = "phate",         # This corresponds to X_phate
  color = "Diagnosis",     # Or "clusters", "core", etc.
  show = FALSE             # Don't display in R; we'll save manually
)

# 3. Save the plot manually using matplotlib
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Diagnosis.png", dpi = 300)



sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Sst_25",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Sst25.png", dpi = 300)



sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Age",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Age.png", dpi = 300)

sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Astro_2",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Astro2.png", dpi = 300)




sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Astro_3",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Astro3.png", dpi = 300)




sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Astro_4",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Astro4.png", dpi = 300)




sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Micro-PVM_2",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Micro2.png", dpi = 300)


sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Micro-PVM_1",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Micro1.png", dpi = 300)



sc$pl$embedding(
   mt_sinai_data,
  basis = "phate",          # assumes you ran PHATE with n_components = 2
  color = "Sst_2",          # or any other column in var_names
  cmap = "viridis",         # optional colormap
  show = FALSE
)

# Save manually
plt <- import("matplotlib.pyplot")
plt$savefig("Figures/Mtsinai_PHATE2D_Sst2.png", dpi = 300)
