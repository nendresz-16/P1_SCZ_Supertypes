setwd ("P1_EDA")
library(Seurat)
library(dplyr)
library(tidyr)
library(readxl)
library(tibble)  
library(dplyr)
library(ggplot2)


#Determine if disease is biased for certain supertypes 
HBCC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/HBCC_updated.rds")
MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
OFC   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds")
Bat   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds")
Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds")

HBCC_meta  <- HBCC@meta.data
MSSM_meta <- rbind(
  MSSM1@meta.data,
  MSSM2@meta.data
)
OFC_meta   <- OFC@meta.data
Bat_meta   <- Bat@meta.data
Ruz_meta   <- Ruz@meta.data




# Extract only the relevant columns from each dataset
HBCC_subset <- HBCC_meta %>% select(predicted.id, Diagnosis, cohort)
MSSM_subset <- MSSM_meta %>% select(predicted.id, Diagnosis, cohort)
OFC_subset  <- OFC_meta  %>% select(predicted.id, Diagnosis, cohort)
Bat_subset  <- Bat_meta  %>% select(predicted.id, Diagnosis, cohort)
Ruz_subset  <- Ruz_meta  %>% select(predicted.id, Diagnosis, cohort)

# Combine all subsets into one dataframe
all_meta <- bind_rows(HBCC_subset, MSSM_subset, OFC_subset, Bat_subset, Ruz_subset)



# Filter for Sst_ predicted cells only
sst_meta <- all_meta %>%
  filter(grepl("^Sst_", predicted.id))

# Count cells per Cohort, Diagnosis, and Sst subtype
sst_counts <- sst_meta %>%
  group_by(cohort, Diagnosis, predicted.id) %>%
  summarise(count = n(), .groups = "drop")

# Plot: one facet per Sst subtype
p <- ggplot(sst_counts, aes(x = cohort, y = count, fill = Diagnosis)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ predicted.id, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Sst Subtypes by Diagnosis Across Cohorts",
    x = "Cohort",
    y = "Cell Count"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold")
  )


ggsave("Figures/Sst_cells_by_Diagnosis.png", p, width = 10, height = 6, dpi=300)