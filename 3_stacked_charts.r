
setwd("P1_Compositional_analysis")
library(Seurat)
library(dplyr)
library(ggplot2)

HBCC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/HBCC_updated.rds")
MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
# Merge

MSSM <- merge(MSSM1, y = MSSM2)
MSSM <- JoinLayers(MSSM)
# ==== 1. Load and subset cohorts ====
OFC   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") 
Bat   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") 
Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 

Ruz_McLean <- subset(Ruz, subset = cohort == "McLean")
Ruz_MtSinai   <- subset(Ruz, subset = cohort == "MtSinai")


# ---- Helper function ----
make_prop_table <- function(seu, cohort_name) {
  seu@meta.data %>%
    count(Donor, predicted.id) %>%
    group_by(Donor) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup() %>%
    mutate(cohort = cohort_name)
}

# ---- Make per-cohort data frames ----
df_HBCC   <- make_prop_table(HBCC, "HBCC")
df_MSSM   <- make_prop_table(MSSM, "MSSM")
df_OFC    <- make_prop_table(OFC, "OFC")
df_Bat    <- make_prop_table(Bat, "Bat")
df_McLean <- make_prop_table(Ruz_McLean, "McLean")
df_MtSinai <- make_prop_table(Ruz_MtSinai, "MtSinai")

# ---- Build data for each cohort ----
all_df <- bind_rows(
  df_HBCC,
  df_MSSM,
  df_OFC,
  df_Bat,
  df_McLean,
  df_MtSinai
)



# ---- Plot ----
# Make plot object
p <- ggplot(all_df, aes(x = Donor, y = freq, fill = predicted.id)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cohort, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Donor",
    y = "Proportion of cells",
    fill = "predicted.id",
    title = "Stacked bar charts of predicted.id per donor across cohorts"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save to file
ggsave("Figures/predicted_id_stacked_bars.png", p,
       width = 30, height = 20, dpi = 300)


###### Percentages

all_df <- all_df %>%
  mutate(predicted.id_family = sub("_.*", "", predicted.id))


# Collapse across donors to get cohort-level percentages
cohort_summary <- all_df %>%
  group_by(cohort, predicted.id_family) %>%
  summarise(total_cells = sum(n), .groups = "drop") %>%
  group_by(cohort) %>%
  mutate(percentage = total_cells / sum(total_cells) * 100)


p1 <- ggplot(cohort_summary, aes(x = "", y = percentage, fill = predicted.id_family)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ cohort) +
  theme_classic() +
  labs(fill = "predicted.id_family", title = "Subclass composition per cohort")

ggsave("Figures/predicted_id_percentage.png", p1,width = 30, height = 20, dpi = 300)





cohort_summary_no_oligo <- all_df %>%
  filter(!(cohort %in% c("HBCC", "MSSM") & grepl("^Oligo", predicted.id_family))) %>%
  group_by(cohort, predicted.id_family) %>%
  summarise(total_cells = sum(n), .groups = "drop") %>%
  group_by(cohort) %>%
  mutate(percentage = total_cells / sum(total_cells) * 100)

# Replot without Oligos
p2 <- ggplot(cohort_summary_no_oligo, aes(x = "", y = percentage, fill = predicted.id_family)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ cohort) +
  theme_classic() +
  labs(fill = "predicted.id_family", title = "Subclass composition per cohort (no Oligos in HBCC/MSSM)")

ggsave("Figures/predicted_id_percentage_noOligo.png", p2, width = 30, height = 20, dpi = 300)






# Keep only SST supertype
sst_summary <- all_df %>%
  filter(grepl("SST", predicted.id, ignore.case = TRUE)) %>%
  group_by(cohort, predicted.id) %>%
  summarise(total_cells = sum(n), .groups = "drop") %>%
  group_by(cohort) %>%
  mutate(percentage = total_cells / sum(total_cells) * 100)

# Plot
p1 <- ggplot(sst_summary, aes(x = "", y = percentage, fill = predicted.id)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ cohort) +
  theme_classic() +
  labs(fill = "predicted.id", title = "SST composition per cohort")

ggsave("Figures/predicted_id_percentage_SST.png", p1,
       width = 30, height = 20, dpi = 300)




family_summary <- all_df %>%
  group_by(cohort, predicted.id_family) %>%
  summarise(total_cells = sum(n), .groups = "drop")




donor_props_family <- all_df %>%
  mutate(predicted.id_family = sub("_.*", "", predicted.id)) %>%
  group_by(cohort, Donor, predicted.id_family) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(cohort, Donor) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

library(broom)

library(ggpubr)

cohorts <- unique(donor_props_family$cohort)

# generate all pairwise combinations
all_pairs <- combn(cohorts, 2, simplify = FALSE)

p3 <- ggplot(donor_props_family, aes(x = cohort, y = prop, fill = cohort)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~ predicted.id_family, scales = "free_y") +
  theme_classic() +
  labs(
    y = "Proportion of cells",
    x = "Cohort",
    fill = "Cohort",
    title = "Predicted.id family proportions across cohorts"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "t.test",
                     comparisons = all_pairs,
                     label = "p.signif")



ggsave("Figures/subclass_ttest.png", p3, width = 30, height = 20, dpi = 300)
