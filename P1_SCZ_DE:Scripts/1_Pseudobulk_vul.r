#conda activate r_env

library(Seurat) 
library(dplyr)

setwd("P1_SCZ_DE")


CELL_TYPES <- c("Sst_1","Sst_4","Sst_5","Sst_7","Sst_3",
                "Sst_19","Sst_9","Sst_10","Sst_13","Sst_12",
                "Sst_11","Sst_20","Sst_22","Sst_23","Sst_25","Sst_2")

affected   <- c("Sst_3","Sst_20","Sst_22","Sst_25","Sst_2", "Sst_23", "Sst_11")
unaffected <- setdiff(CELL_TYPES, affected)


OFC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") 

Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") 

Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 

# Split Ruz into McLean and MtSinai
Ruz_McLean  <- subset(Ruz, subset = Cohort == "McLean")
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


#Cell vulnerability 




# 1) Subset to the SST subtypes of interest
OFC_sst <- subset(OFC, subset = predicted.id %in% CELL_TYPES)

# 2) Add binary label
OFC_sst$Vulnerability <- ifelse(OFC_sst$predicted.id %in% affected, "Affected", "Unaffected")

# 3) Aggregate to pseudobulk (Vulnerability × Donor)
bulk_OFC <- AggregateExpression(OFC_sst,  group.by = c("Vulnerability", "Donor"),  return.seurat = TRUE)

raw_counts_OFC <- GetAssayData(object = bulk_OFC, assay = "RNA", layer = "counts")

saveRDS(raw_counts_OFC, "Files/Pseudobulk_OFC.rds")


# --- Bat ---
Bat_sst <- subset(Bat, subset = predicted.id %in% CELL_TYPES)
Bat_sst$Vulnerability <- ifelse(Bat_sst$predicted.id %in% affected, "Affected", "Unaffected")
bulk_Bat <- AggregateExpression(Bat_sst, group.by = c("Vulnerability", "Donor"), return.seurat = TRUE)
raw_counts_Bat <- GetAssayData(object = bulk_Bat, assay = "RNA", layer = "counts")
saveRDS(raw_counts_Bat, "Files/Pseudobulk_Bat.rds")

# --- Ruz_McLean ---
Ruz_McLean_sst <- subset(Ruz_McLean, subset = predicted.id %in% CELL_TYPES)
Ruz_McLean_sst$Vulnerability <- ifelse(Ruz_McLean_sst$predicted.id %in% affected, "Affected", "Unaffected")
bulk_McLean <- AggregateExpression(Ruz_McLean_sst, group.by = c("Vulnerability", "Donor"), return.seurat = TRUE)
raw_counts_McLean <- GetAssayData(object = bulk_McLean, assay = "RNA", layer = "counts")
saveRDS(raw_counts_McLean, "Files/Pseudobulk_McLean.rds")

# --- Ruz_MtSinai ---
Ruz_MtSinai_sst <- subset(Ruz_MtSinai, subset = predicted.id %in% CELL_TYPES)
Ruz_MtSinai_sst$Vulnerability <- ifelse(Ruz_MtSinai_sst$predicted.id %in% affected, "Affected", "Unaffected")
bulk_MtSinai <- AggregateExpression(Ruz_MtSinai_sst, group.by = c("Vulnerability", "Donor"), return.seurat = TRUE)
raw_counts_MtSinai <- GetAssayData(object = bulk_MtSinai, assay = "RNA", layer = "counts")
saveRDS(raw_counts_MtSinai, "Files/Pseudobulk_MtSinai.rds")






#Combined metadata 
vars_to_keep <- c("Donor", "Diagnosis", "Sex", "cohort", "predicted.id", "Age", "PMI", "Vulnerability")

meta_ofc     <- OFC_sst@meta.data     %>% dplyr::select(any_of(vars_to_keep))
meta_bat     <- Bat_sst@meta.data     %>% dplyr::select(any_of(vars_to_keep))
meta_mclean  <- Ruz_McLean_sst@meta.data %>% dplyr::select(any_of(vars_to_keep))
meta_mtsinai <- Ruz_MtSinai_sst@meta.data %>% dplyr::select(any_of(vars_to_keep))

metadata_all <- bind_rows(meta_ofc, meta_bat, meta_mclean, meta_mtsinai)


#Pseudobulk meta data

# Step 1: Count number of cells per Donor × Supertypes

donor_counts <- metadata_all %>%

group_by(Vulnerability, Donor) %>%

summarise(num_cells = n(), .groups = "drop")

  

# Step 2: Extract donor-level metadata 

donor_metadata <- metadata_all%>%

group_by(Donor) %>% summarise(Age = first(Age), Sex = first(Sex), PMI = first(PMI), Diagnosis = first(Diagnosis), cohort = first(cohort), .groups = "drop")

  

# Step 3: Combine cell counts with donor metadata

pseudobulk_meta <- left_join(donor_metadata, donor_counts, by = "Donor")

write.csv(pseudobulk_meta, "Files/pseudobulk_metadata_SCZ.csv", row.names = FALSE)