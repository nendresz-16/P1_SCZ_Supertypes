setwd("P1_Compositional_analysis")

HBCC <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/HBCC_updated.rds")
MSSM1 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM1_updated.rds")
MSSM2 <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/MSSM2_updated.rds")
# Merge

MSSM <- merge(MSSM1, y = MSSM2)


OFC   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds") 
Bat   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") 
Ruz   <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") 


# Split Ruz into McLean and MtSinai
Ruz_McLean  <- subset(Ruz, subset = Cohort == "McLean")
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")


#Extract meta data 

hbcc_meta       <- HBCC@meta.data
mssm_meta       <- MSSM@meta.data
ofc_meta        <- OFC@meta.data
bat_meta        <- Bat@meta.data
ruz_mclean_meta <- Ruz_McLean@meta.data
ruz_mtsinai_meta<- Ruz_MtSinai@meta.data


library(dplyr)
library(stringr)
library(tidyr)

# Classify into 3 broad classes
classify_cell <- function(x) {
  if (str_detect(x, "L2/3 IT|L4 IT|L5 IT|L6 IT|L6b|L5 ET|L6 CT|L5/6 NP|Car3")) {
    return("Glutamatergic")
  } else if (str_detect(x, "Sst|Pvalb|Vip|Lamp5|Sncg|Pax6|Chandelier")) {
    return("GABAergic")
  } else {
    return("Non-neuronal")
  }
}

# Summarize directly to donor × broad class
summarize_dataset_broad <- function(meta, dataset_name) {
  meta %>%
    mutate(cell_class = sapply(predicted.id, classify_cell)) %>%
    group_by(Donor, cell_class) %>%
    summarise(n_class = n(), .groups = "drop_last") %>%
    group_by(Donor) %>%
    mutate(
      total_cells = sum(n_class),
      prop_class = n_class / total_cells,
      dataset = dataset_name
    ) %>%
    ungroup() %>%
    select(Donor, dataset, total_cells, cell_class, prop_class) %>%
    pivot_wider(
      names_from = cell_class,
      values_from = prop_class,
      values_fill = 0
    )
}

# Run for all cohorts
hbcc_summary_broad       <- summarize_dataset_broad(hbcc_meta, "HBCC")
mssm_summary_broad       <- summarize_dataset_broad(mssm_meta, "MSSM")
ofc_summary_broad        <- summarize_dataset_broad(ofc_meta, "OFC")
bat_summary_broad        <- summarize_dataset_broad(bat_meta, "Bat")
ruz_mclean_summary_broad <- summarize_dataset_broad(ruz_mclean_meta, "Ruz_McLean")
ruz_mtsinai_summary_broad<- summarize_dataset_broad(ruz_mtsinai_meta, "Ruz_MtSinai")



QC_MIN_CELLS   <- 1500   # or 2000 if you want Mathys criteria
QC_NN_PROP_MIN <- 0.125
QC_NN_PROP_MAX <- 0.67


ruz_mtsinai_qc <- ruz_mtsinai_summary_broad %>%
  rename(prop_non_neuronal = `Non-neuronal`) %>%
  mutate(
    qc_status = case_when(
      total_cells < QC_MIN_CELLS ~ "failed_total_cells",
      prop_non_neuronal < QC_NN_PROP_MIN ~ "failed_nn_prop_low",
      prop_non_neuronal > QC_NN_PROP_MAX ~ "failed_nn_prop_high",
      TRUE ~ "PASS"
    )
  )


library(ggplot2)

p_ruz <- ggplot(ruz_mtsinai_qc, aes(x = total_cells, y = prop_non_neuronal, color = qc_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = QC_MIN_CELLS, linetype = "dashed", color = "brown") +
  geom_hline(yintercept = QC_NN_PROP_MIN, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = QC_NN_PROP_MAX, linetype = "dotted", color = "blue") +
  theme_classic() +
  labs(
    title = "Ruz_MtSinai QC Verification",
    x = "Total Cells per Donor",
    y = "Proportion of Non-Neuronal Cells"
  ) +
  scale_color_manual(
    values = c(
      "failed_nn_prop_high" = "#F8766D",   # reddish
      "failed_total_cells" = "#00BA38",    # green
      "failed_nn_prop_low" = "#B79F00",    # yellow-brown
      "failed_total_cells; failed_nn_prop_low" = "#619CFF", # blue
      "failed_total_cells; failed_nn_prop_high" = "#00BFC4", # teal
      "PASS" = "#E76BF3"                   # magenta/pink
    )
  )

# Save as PNG
ggsave("Figures/Ruz_MtSinai_QC.png", p_ruz, width = 6, height = 5, dpi = 300)



hbcc_qc <- hbcc_summary_broad %>%
  rename(prop_non_neuronal = `Non-neuronal`) %>%
  mutate(
    qc_status = case_when(
      total_cells < QC_MIN_CELLS ~ "failed_total_cells",
      prop_non_neuronal < QC_NN_PROP_MIN ~ "failed_nn_prop_low",
      prop_non_neuronal > QC_NN_PROP_MAX ~ "failed_nn_prop_high",
      TRUE ~ "PASS"
    )
  )

p_hbcc <- ggplot(hbcc_qc, aes(x = total_cells, y = prop_non_neuronal, color = qc_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = QC_MIN_CELLS, linetype = "dashed", color = "brown") +
  geom_hline(yintercept = QC_NN_PROP_MIN, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = QC_NN_PROP_MAX, linetype = "dotted", color = "blue") +
  theme_classic() +
  labs(
    title = "HBCC QC Verification",
    x = "Total Cells per Donor",
    y = "Proportion of Non-Neuronal Cells"
  )+
  scale_color_manual(
    values = c(
      "failed_nn_prop_high" = "#F8766D",   # reddish
      "failed_total_cells" = "#00BA38",    # green
      "failed_nn_prop_low" = "#B79F00",    # yellow-brown
      "failed_total_cells; failed_nn_prop_low" = "#619CFF", # blue
      "failed_total_cells; failed_nn_prop_high" = "#00BFC4", # teal
      "PASS" = "#E76BF3"                   # magenta/pink
    )
  )

ggsave("Figures/HBCC_QC.png", p_hbcc, width = 6, height = 5, dpi = 300)

mssm_qc <- mssm_summary_broad %>%
  rename(prop_non_neuronal = `Non-neuronal`) %>%
  mutate(
    qc_status = case_when(
      total_cells < QC_MIN_CELLS ~ "failed_total_cells",
      prop_non_neuronal < QC_NN_PROP_MIN ~ "failed_nn_prop_low",
      prop_non_neuronal > QC_NN_PROP_MAX ~ "failed_nn_prop_high",
      TRUE ~ "PASS"
    )
  )

p_mssm <- ggplot(mssm_qc, aes(x = total_cells, y = prop_non_neuronal, color = qc_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = QC_MIN_CELLS, linetype = "dashed", color = "brown") +
  geom_hline(yintercept = QC_NN_PROP_MIN, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = QC_NN_PROP_MAX, linetype = "dotted", color = "blue") +
  theme_classic() +
  labs(
    title = "MSSM QC Verification",
    x = "Total Cells per Donor",
    y = "Proportion of Non-Neuronal Cells"
  ) +
  scale_color_manual(
    values = c(
      "failed_nn_prop_high" = "#F8766D",   # reddish
      "failed_total_cells" = "#00BA38",    # green
      "failed_nn_prop_low" = "#B79F00",    # yellow-brown
      "failed_total_cells; failed_nn_prop_low" = "#619CFF", # blue
      "failed_total_cells; failed_nn_prop_high" = "#00BFC4", # teal
      "PASS" = "#E76BF3"                   # magenta/pink
    )
  )

ggsave("Figures/MSSM_QC.png", p_mssm, width = 6, height = 5, dpi = 300)


ofc_qc <- ofc_summary_broad %>%
  rename(prop_non_neuronal = `Non-neuronal`) %>%
  mutate(
    qc_status = case_when(
      total_cells < QC_MIN_CELLS ~ "failed_total_cells",
      prop_non_neuronal < QC_NN_PROP_MIN ~ "failed_nn_prop_low",
      prop_non_neuronal > QC_NN_PROP_MAX ~ "failed_nn_prop_high",
      TRUE ~ "PASS"
    )
  )

p_ofc <- ggplot(ofc_qc, aes(x = total_cells, y = prop_non_neuronal, color = qc_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = QC_MIN_CELLS, linetype = "dashed", color = "brown") +
  geom_hline(yintercept = QC_NN_PROP_MIN, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = QC_NN_PROP_MAX, linetype = "dotted", color = "blue") +
  theme_classic() +
  labs(
    title = "OFC QC Verification",
    x = "Total Cells per Donor",
    y = "Proportion of Non-Neuronal Cells"
  ) +
  scale_color_manual(
    values = c(
      "failed_nn_prop_high" = "#F8766D",   # reddish
      "failed_total_cells" = "#00BA38",    # green
      "failed_nn_prop_low" = "#B79F00",    # yellow-brown
      "failed_total_cells; failed_nn_prop_low" = "#619CFF", # blue
      "failed_total_cells; failed_nn_prop_high" = "#00BFC4", # teal
      "PASS" = "#E76BF3"                   # magenta/pink
    )
  )

ggsave("Figures/OFC_QC.png", p_ofc, width = 6, height = 5, dpi = 300)



bat_qc <- bat_summary_broad %>%
  rename(prop_non_neuronal = `Non-neuronal`) %>%
  mutate(
    qc_status = case_when(
      total_cells < QC_MIN_CELLS ~ "failed_total_cells",
      prop_non_neuronal < QC_NN_PROP_MIN ~ "failed_nn_prop_low",
      prop_non_neuronal > QC_NN_PROP_MAX ~ "failed_nn_prop_high",
      TRUE ~ "PASS"
    )
  )

library(ggplot2)

p_bat <- ggplot(bat_qc, aes(x = total_cells, y = prop_non_neuronal, color = qc_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = QC_MIN_CELLS, linetype = "dashed", color = "brown") +
  geom_hline(yintercept = QC_NN_PROP_MIN, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = QC_NN_PROP_MAX, linetype = "dotted", color = "blue") +
  theme_classic() +
  labs(
    title = "Bat QC Verification",
    x = "Total Cells per Donor",
    y = "Proportion of Non-Neuronal Cells"
  ) +
  scale_color_manual(
    values = c(
      "failed_nn_prop_high" = "#F8766D",   # reddish
      "failed_total_cells" = "#00BA38",    # green
      "failed_nn_prop_low" = "#B79F00",    # yellow-brown
      "failed_total_cells; failed_nn_prop_low" = "#619CFF", # blue
      "failed_total_cells; failed_nn_prop_high" = "#00BFC4", # teal
      "PASS" = "#E76BF3"                   # magenta/pink
    )
  )

# Save as PNG
ggsave("Figures/Bat_QC.png", p_bat, width = 6, height = 5, dpi = 300)



ruz_mclean_qc <- ruz_mclean_summary_broad %>%
  rename(prop_non_neuronal = `Non-neuronal`) %>%
  mutate(
    qc_status = case_when(
      total_cells < QC_MIN_CELLS ~ "failed_total_cells",
      prop_non_neuronal < QC_NN_PROP_MIN ~ "failed_nn_prop_low",
      prop_non_neuronal > QC_NN_PROP_MAX ~ "failed_nn_prop_high",
      TRUE ~ "PASS"
    )
  )

p_ruz_mclean <- ggplot(ruz_mclean_qc, aes(x = total_cells, y = prop_non_neuronal, color = qc_status)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = QC_MIN_CELLS, linetype = "dashed", color = "brown") +
  geom_hline(yintercept = QC_NN_PROP_MIN, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = QC_NN_PROP_MAX, linetype = "dotted", color = "blue") +
  theme_classic() +
  labs(
    title = "Ruz_McLean QC Verification",
    x = "Total Cells per Donor",
    y = "Proportion of Non-Neuronal Cells"
  ) +
  scale_color_manual(
    values = c(
      "failed_nn_prop_high" = "#F8766D",   # reddish
      "failed_total_cells" = "#00BA38",    # green
      "failed_nn_prop_low" = "#B79F00",    # yellow-brown
      "failed_total_cells; failed_nn_prop_low" = "#619CFF", # blue
      "failed_total_cells; failed_nn_prop_high" = "#00BFC4", # teal
      "PASS" = "#E76BF3"                   # magenta/pink
    )
  )

# Save as PNG
ggsave("Figures/Ruz_McLean_QC.png", p_ruz_mclean, width = 6, height = 5, dpi = 300)



library(patchwork)

# Arrange in a 2 x 3 grid
combined <- (p_hbcc + p_mssm + p_ofc) /
            (p_bat + p_ruz_mclean + p_ruz)

# Show in R
combined

# Save combined figure
ggsave("Figures/All_Cohorts_QC.png", combined, width = 18, height = 12, dpi = 300)
