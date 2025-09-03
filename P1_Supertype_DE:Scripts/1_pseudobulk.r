#conda activate r_env

library(Seurat) 
library(dplyr)
setwd("P1_Supertype_DE")

obj <- readRDS("/project/rrg-shreejoy/nendresz/SEAD_data/SEAD_obj_gene_symbols.rds")

levels(factor(obj@meta.data$Cognitive.status))


#Filter out individuals w dementia 

obj <- subset(obj, subset = Cognitive.status %in% c("Reference", "No dementia"))

#How many donors 

length(unique(obj@meta.data$donor_id))

#Pseudobulk counts by supertype
bulk <- AggregateExpression(obj, group.by = c("Supertype", "donor_id"), return.seurat = TRUE)

  
Cells(bulk)

  

head(bulk$RNA$counts) # raw

  

raw_counts <- GetAssayData(object = bulk, assay = "RNA", layer = "counts")


saveRDS(raw_counts, "Files/Pseudobulk_SEAD_con.rds")


#Pseudobulk meta data

# Step 1: Count number of cells per Donor × Supertypes

donor_counts <- obj@meta.data %>%

group_by(donor_id, Supertype) %>%

summarise(num_cells = n(), .groups = "drop")

  

# Step 2: Extract donor-level metadata 

donor_metadata <- obj@meta.data %>%

group_by(donor_id) %>% summarise(Age = first(Age.at.death), Sex = first(sex), PMI = first(PMI), .groups = "drop")

  

# Step 3: Combine cell counts with donor metadata

pseudobulk_SEAD <- left_join(donor_metadata, donor_counts, by = "donor_id")

write.csv(pseudobulk_SEAD, "Files/pseudobulk_metadata_SEAD_con.csv", row.names = FALSE)