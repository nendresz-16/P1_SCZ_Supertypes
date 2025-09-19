setwd("P1_Trajectory_analysis")

library(Seurat)

Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") 
OFC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") 
Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 

# Split Ruz into McLean and MtSinai
Ruz_McLean  <- subset(Ruz, subset = Cohort == "McLean")
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")



###################
###### BATIUK #####
###################
# genes x cells sparse matrix
counts <- GetAssayData(Bat, assay = "RNA", layer = "counts")

# Write Matrix Market + gene/cell lists
Matrix::writeMM(counts, "Files/Bat_counts.mtx")
writeLines(rownames(counts), "Files/Bat_genes.txt")
writeLines(colnames(counts), "Files/Bat_cells.txt")

# Metadata (cell-level)
meta <- Bat@meta.data
meta$cell <- rownames(meta)
write.csv(meta, "Files/Bat_metadata.csv", row.names = FALSE)


