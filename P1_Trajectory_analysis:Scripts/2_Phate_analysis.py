# conda activate phate_env


import pandas as pd, numpy as np
import anndata as ad
import os
import scanpy as sc
import phate



base = "/scratch/nendresz/P1_Trajectory_analysis/Files"

# Load sparse counts (genes x cells), then transpose -> cells x genes
X = sio.mmread(f"{base}/Bat_counts.mtx").tocsr().T
genes = pd.read_csv(f"{base}/Bat_genes.txt", header=None, sep="\t").iloc[:,0].astype(str).values
cells = pd.read_csv(f"{base}/Bat_cells.txt", header=None, sep="\t").iloc[:,0].astype(str).values

meta = pd.read_csv(f"{base}/Bat_metadata.csv", index_col="cell", dtype=str)
meta = meta.loc[cells]  # keep same order as matrix

adata = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
adata.obs_names_make_unique(); adata.var_names_make_unique()
adata.write_h5ad(f"{base}/Bat.h5ad", compression="lzf")

####
import anndata as ad

adata = ad.read_h5ad("/scratch/nendresz/P1_Trajectory_analysis/Files/Bat.h5ad")

import phate
phate_op = phate.PHATE()
data_phate = phate_op.fit_transform(adata)
######

import pandas as pd
import numpy as np
import phate

# 0) pick the donor id column and a cell-type label column from adata.obs
don_col = "Donor"           # or "DonorLabel"
ct_col  = "predicted.id"    # or your preferred cell-type/cluster label

# 1) donor × cell-type proportions
tab  = pd.crosstab(adata.obs[don_col], adata.obs[ct_col])
prop = tab.div(tab.sum(axis=1), axis=0).fillna(0)   # normalize rows to proportions

# ---------- ALL DONORS: 3D PHATE ----------
ph_all = phate.PHATE(n_components=3, k=10, a=40, knn_dist='euclidean', mds_dist='euclidean', n_jobs=-1)
Y_all  = ph_all.fit_transform(prop.values)
df_all = pd.DataFrame(Y_all, index=prop.index, columns=["PHATE1","PHATE2","PHATE3"])

# store somewhere handy
# aligned to donors (not cells):
donor_phate_all = df_all

# ---------- CORE DONORS: 2D PHATE ----------
# assumes a boolean flag in adata.obs called 'core' (True/False)
core_donors = adata.obs.loc[adata.obs["core"].astype(bool), don_col].unique()
prop_core   = prop.loc[core_donors]

ph_core = phate.PHATE(n_components=2, k=15, a=100, knn_dist='euclidean', mds_dist='correlation', n_jobs=-1)
Y_core  = ph_core.fit_transform(prop_core.values)
df_core = pd.DataFrame(Y_core, index=prop_core.index, columns=["PHATE1","PHATE2"])

# make an all-donor–indexed frame with NaNs for non-core donors (like their merge)
df_core_aligned = pd.DataFrame(index=prop.index, columns=df_core.columns)
df_core_aligned.loc[df_core.index] = df_core.values

# OPTIONAL: broadcast donor embeddings back to cells so you can color cells by donor-PHATE
# (creates cell-level arrays aligned to adata.obs_names)
adata.obsm["X_donor_all_3d_phate"] = donor_phate_all.reindex(adata.obs[don_col]).values
adata.obsm["X_donor_core_2d_phate"] = df_core_aligned.reindex(adata.obs[don_col]).values
