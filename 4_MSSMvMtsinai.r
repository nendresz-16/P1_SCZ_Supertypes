setwd("P1_Compositional_analysis")
library(Seurat)
library(dplyr)
library(ggplot2)
library(MatchIt)


MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
MSSM <- merge(MSSM1, y = MSSM2)
MSSM <- JoinLayers(MSSM)

Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


# Collapse MSSM donors
MS_donors <- MSSM@meta.data %>%
  select(Donor, Sex, Diagnosis, Age) %>%
  distinct() %>%
  mutate(Cohort = "MSSM")

# Collapse MtSinai donors
MT_donors <- Ruz_MtSinai@meta.data %>%
  select(Donor, Sex, Diagnosis, Age) %>%
  distinct() %>%
  mutate(Cohort = "MtSinai")

# Combine donor-level metadata
combined <- bind_rows(MS_donors, MT_donors)


library(MatchIt)

combined$Cohort <- factor(combined$Cohort, levels = c("MSSM", "MtSinai"))

m.out <- matchit(
  Cohort ~ Age + Sex + Diagnosis,
  data = combined,
  method = "nearest"
)

matched_data <- match.data(m.out)




pairs <- matched_data %>%
  group_by(subclass) %>%
  summarise(
    MSSM_donor       = Donor[Cohort == "MSSM"],
    MtSinai_donor    = Donor[Cohort == "MtSinai"],
    Age_MSSM         = Age[Cohort == "MSSM"],
    Age_MtSinai      = Age[Cohort == "MtSinai"],
    Sex_MSSM         = Sex[Cohort == "MSSM"],
    Sex_MtSinai      = Sex[Cohort == "MtSinai"],
    Diagnosis_MSSM   = Diagnosis[Cohort == "MSSM"],
    Diagnosis_MtSinai= Diagnosis[Cohort == "MtSinai"],
    .groups = "drop"
  )

head(pairs)

write.csv(pairs, "Files/MSSM_MtSinai_matched_pairs.csv", row.names = FALSE)





# Find exact matches on Age, Sex, and Diagnosis
perfect_matches <- pairs %>%
  filter(
    Age_MSSM == Age_MtSinai,
    Sex_MSSM == Sex_MtSinai,
    Diagnosis_MSSM == Diagnosis_MtSinai
  )

# Show just the donors if you want
perfect_donors <- perfect_matches %>%
  select(MSSM_donor, MtSinai_donor)

print(perfect_donors)

perfect_donors <- perfect_donors %>%
  mutate(pair_id = paste0("Pair_", row_number()))

print(perfect_donors)





#####################################

MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
MSSM <- merge(MSSM1, y = MSSM2)
MSSM <- JoinLayers(MSSM)

Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


# Take the part before the first "_" in predicted.id
MSSM$subclass_type        <- sub("_.*", "", MSSM$predicted.id)
Ruz_MtSinai$subclass_type <- sub("_.*", "", Ruz_MtSinai$predicted.id)


MSSM_prop <- MSSM@meta.data %>%
  count(cohort, Donor, subclass_type) %>%
  group_by(cohort, Donor) %>%
  mutate(total_cells = sum(n)) %>%
  ungroup() %>%
  filter(total_cells >= 500) %>%        # drop donors with <500 cells
  group_by(cohort, Donor) %>%
  mutate(prop = n / total_cells) %>%
  ungroup()

  Ruz_props <- Ruz_MtSinai@meta.data %>%
  count(cohort, Donor, subclass_type) %>%
  group_by(cohort, Donor) %>%
  mutate(total_cells = sum(n)) %>%
  ungroup() %>%
  filter(total_cells >= 500) %>%        # drop donors with <500 cells
  group_by(cohort, Donor) %>%
  mutate(prop = n / total_cells) %>%
  ungroup()

  # --- 2. Convert to wide (per donor × subclass_type) ---
df_mssm <- MSSM_prop %>%
  select(Donor, subclass_type, prop) %>%
  pivot_wider(names_from = subclass_type, values_from = prop, values_fill = 0)

df_mt <- Ruz_props %>%
  select(Donor, subclass_type, prop) %>%
  pivot_wider(names_from = subclass_type, values_from = prop, values_fill = 0)

# --- 3. Merge matched pairs ---
pairs_props <- perfect_donors %>%
  left_join(df_mssm,  by = c("MSSM_donor"   = "Donor")) %>%
  left_join(df_mt,    by = c("MtSinai_donor" = "Donor"),
            suffix = c("_MSSM", "_MtSinai"))

# --- 4. Reshape for plotting ---
long_pairs <- pairs_props %>%
  pivot_longer(
    cols = matches("_(MSSM|MtSinai)$"),
    names_to = c("subclass_type", "cohort"),
    names_pattern = "(.*)_(MSSM|MtSinai)",
    values_to = "prop"
  ) %>%
  pivot_wider(names_from = cohort, values_from = prop)

  # --- Filter out L5 ET and Sst Chodl before plotting ---
long_pairs_filtered <- long_pairs %>%
  filter(!subclass_type %in% c("L5 ET", "Sst Chodl"))

# --- 5. Scatterplot ---
p_all <- ggplot(long_pairs_filtered, aes(x = MSSM, y = MtSinai, label = pair_id)) +
  geom_point(size = 2, alpha = 0.8, color = "steelblue") +
  geom_text(vjust = -0.5, size = 2.5) +   # annotate with pair IDs
  geom_smooth(method = "lm", se = FALSE, color = "blue") + # regression line only
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  facet_wrap(~ subclass_type, scales = "free") +
  theme_classic() +
  labs(
    x = "MSSM proportion",
    y = "MtSinai proportion",
    title = "Matched donor subclass proportions (MSSM vs MtSinai)"
  )

ggsave("Figures/MSSM_vs_MtSinai_matched_pairs_filtered.png", p_all,
       width = 16, height = 10, dpi = 300)


################################

MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
MSSM <- merge(MSSM1, y = MSSM2)
MSSM <- JoinLayers(MSSM)

Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


# === 3. Add subclass_type (prefix of predicted.id) ===
MSSM$subclass_type        <- sub("_.*", "", MSSM$predicted.id)
Ruz_MtSinai$subclass_type <- sub("_.*", "", Ruz_MtSinai$predicted.id)


# === 4. Keep neurons only ===
neuron_types <- c(
  "L2/3 IT", "L4 IT", "L5 IT", "L5/6 NP",
  "L6 IT", "L6 IT Car3", "L6 CT", "L6b",
  "Sst", "Pvalb", "Vip", "Lamp5", "Sncg", "Chandelier", "L5 ET", "Sst Chodl"
)

MSSM  <- subset(MSSM, subset = subclass_type %in% neuron_types)
Ruz_MtSinai <- subset(Ruz_MtSinai, subset = subclass_type %in% neuron_types)


MSSM_prop <- MSSM@meta.data %>%
  count(cohort, Donor, subclass_type) %>%
  group_by(cohort, Donor) %>%
  mutate(total_cells = sum(n)) %>%
  ungroup() %>%
  group_by(cohort, Donor) %>%
  mutate(prop = n / total_cells) %>%
  ungroup()

  Ruz_props <- Ruz_MtSinai@meta.data %>%
  count(cohort, Donor, subclass_type) %>%
  group_by(cohort, Donor) %>%
  mutate(total_cells = sum(n)) %>%
  ungroup() %>%
  group_by(cohort, Donor) %>%
  mutate(prop = n / total_cells) %>%
  ungroup()

  # --- 2. Convert to wide (per donor × subclass_type) ---
df_mssm <- MSSM_prop %>%
  select(Donor, subclass_type, prop) %>%
  pivot_wider(names_from = subclass_type, values_from = prop, values_fill = 0)

df_mt <- Ruz_props %>%
  select(Donor, subclass_type, prop) %>%
  pivot_wider(names_from = subclass_type, values_from = prop, values_fill = 0)

# --- 3. Merge matched pairs ---
pairs_props <- perfect_donors %>%
  left_join(df_mssm,  by = c("MSSM_donor"   = "Donor")) %>%
  left_join(df_mt,    by = c("MtSinai_donor" = "Donor"),
            suffix = c("_MSSM", "_MtSinai"))

# --- 4. Reshape for plotting ---
long_pairs <- pairs_props %>%
  pivot_longer(
    cols = matches("_(MSSM|MtSinai)$"),
    names_to = c("subclass_type", "cohort"),
    names_pattern = "(.*)_(MSSM|MtSinai)",
    values_to = "prop"
  ) %>%
  pivot_wider(names_from = cohort, values_from = prop)

long_pairs_filtered <- long_pairs %>%
  filter(!subclass_type %in% c("L5 ET", "Sst Chodl"))

p_all <- ggplot(long_pairs_filtered, aes(x = MSSM, y = MtSinai, label = pair_id)) +
  geom_point(size = 2, alpha = 0.8, color = "steelblue") +
  geom_text(vjust = -0.5, size = 2.5) +   # annotate with pair IDs
  geom_smooth(method = "lm", se = FALSE, color = "blue") + # regression line only
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  facet_wrap(~ subclass_type, scales = "free") +
  theme_classic() +
  labs(
    x = "MSSM proportion",
    y = "MtSinai proportion",
    title = "Matched donor subclass proportions (MSSM vs MtSinai)"
  )

ggsave("Figures/MSSM_vs_MtSinai_Neurons.png", p_all,
       width = 16, height = 10, dpi = 300)



       ######################


MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
MSSM <- merge(MSSM1, y = MSSM2)
MSSM <- JoinLayers(MSSM)

Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


# === 3. Add subclass_type (prefix of predicted.id) ===
MSSM$subclass_type        <- sub("_.*", "", MSSM$predicted.id)
Ruz_MtSinai$subclass_type <- sub("_.*", "", Ruz_MtSinai$predicted.id)


# === 4. Keep neurons only ===
neuron_types <- c("Sst")

MSSM  <- subset(MSSM, subset = subclass_type %in% neuron_types)
Ruz_MtSinai <- subset(Ruz_MtSinai, subset = subclass_type %in% neuron_types)


MSSM_prop <- MSSM@meta.data %>%
  count(cohort, Donor, predicted.id) %>%
  group_by(cohort, Donor) %>%
  mutate(total_cells = sum(n)) %>%
  ungroup() %>%
  group_by(cohort, Donor) %>%
  mutate(prop = n / total_cells) %>%
  ungroup()

  Ruz_props <- Ruz_MtSinai@meta.data %>%
  count(cohort, Donor, predicted.id) %>%
  group_by(cohort, Donor) %>%
  mutate(total_cells = sum(n)) %>%
  ungroup() %>%
  group_by(cohort, Donor) %>%
  mutate(prop = n / total_cells) %>%
  ungroup()

  # --- 2. Convert to wide (per donor × predicted.id) ---
df_mssm <- MSSM_prop %>%
  select(Donor, predicted.id, prop) %>%
  pivot_wider(names_from = predicted.id, values_from = prop, values_fill = 0)

df_mt <- Ruz_props %>%
  select(Donor, predicted.id, prop) %>%
  pivot_wider(names_from = predicted.id, values_from = prop, values_fill = 0)

# --- 3. Merge matched pairs ---
pairs_props <- perfect_donors %>%
  left_join(df_mssm,  by = c("MSSM_donor"   = "Donor")) %>%
  left_join(df_mt,    by = c("MtSinai_donor" = "Donor"),
            suffix = c("_MSSM", "_MtSinai"))

# --- 4. Reshape for plotting ---
long_pairs <- pairs_props %>%
  pivot_longer(
    cols = matches("_(MSSM|MtSinai)$"),
    names_to = c("predicted.id", "cohort"),
    names_pattern = "(.*)_(MSSM|MtSinai)",
    values_to = "prop"
  ) %>%
  pivot_wider(names_from = cohort, values_from = prop)


p_all <- ggplot(long_pairs, aes(x = MSSM, y = MtSinai, label = pair_id)) +
  geom_point(size = 2, alpha = 0.8, color = "steelblue") +
  geom_text(vjust = -0.5, size = 2.5) +   # annotate with pair IDs
  geom_smooth(method = "lm", se = FALSE, color = "blue") + # regression line only
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  facet_wrap(~ predicted.id, scales = "free") +
  theme_classic() +
  labs(
    x = "MSSM proportion",
    y = "MtSinai proportion",
    title = "Matched donor subclass proportions (MSSM vs MtSinai)"
  )

ggsave("Figures/MSSM_vs_MtSinai_Supertypes.png", p_all,
       width = 16, height = 10, dpi = 300)