# conda activate crumblr_env
# Load libraries
library(SingleCellExperiment)
library(zellkonverter)
library(ggplot2)
library(scattermore)
library(ggtree)
library(crumblr)
library(aplot)
library(tidyverse)
library(dreamlet)
library(kableExtra)
library(ggcorrplot)
library(RColorBrewer)
library(DelayedArray)
library(Seurat)
library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(reshape2)
library(ggplot2)
library(ape)
setwd("P1_Compositional_analysis")



# ==== 1. Load and subset cohorts ( Donor ≥ 500 cells) ====

# Helper function to filter by age and donor cell count
filter_seurat_by_cells <- function(seu, min_cells = 500) {
  donor_counts <- table(seu$Donor)
  donors_to_keep <- names(donor_counts)[donor_counts >= min_cells]
  subset(seu, subset = Donor %in% donors_to_keep)
}

Bat <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds") %>%
  filter_seurat_by_cells()


Bat@meta.data <- Bat@meta.data %>%
  mutate(
    neuroleptic = case_when(
      Typical_neuroleptics == "Y" & Atypical_neuroleptics == "Y" ~ "both",
      Typical_neuroleptics == "Y" & Atypical_neuroleptics == "N" ~ "typ",
      Typical_neuroleptics == "N" & Atypical_neuroleptics == "Y" ~ "atyp",
      Typical_neuroleptics == "N" & Atypical_neuroleptics == "N" ~ "none",
      TRUE ~ NA_character_
    ),
    neuroleptic = factor(neuroleptic, levels = c("none","typ","atyp","both")))

library(dplyr)

donor_neuroleptic <- Bat@meta.data %>%
  distinct(Donor = Donor, neuroleptic,Atypical_neuroleptics,Typical_neuroleptics)  


Ruz <- readRDS("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds") %>%
  filter_seurat_by_cells()

# Split Ruz into McLean and MtSinai
Ruz_McLean  <- subset(Ruz, subset = Cohort == "McLean")
Ruz_MtSinai <- subset(Ruz, subset = Cohort == "MtSinai")



Ruz_McLean@meta.data <- Ruz_McLean@meta.data %>%
  mutate(
    neuroleptic = case_when(
      AntipsychTyp == "1" & AntipsychAtyp == "1" ~ "both",
      AntipsychTyp == "1" & AntipsychAtyp == "0" ~ "typ",
      AntipsychTyp == "0" & AntipsychAtyp== "1" ~ "atyp",
      AntipsychTyp == "0" & AntipsychAtyp == "0" ~ "none",
      TRUE ~ NA_character_
    ),
    neuroleptic = factor(neuroleptic, levels = c("none","typ","atyp","both")))


donor_neuroleptic <- Ruz_McLean@meta.data %>%
  distinct(Donor = Donor, neuroleptic)  



Ruz_MtSinai@meta.data <- Ruz_MtSinai@meta.data %>%
  mutate(
    neuroleptic = case_when(
      AntipsychTyp == "1" & AntipsychAtyp == "1" ~ "both",
      AntipsychTyp == "1" & AntipsychAtyp == "0" ~ "typ",
      AntipsychTyp == "0" & AntipsychAtyp== "1" ~ "atyp",
      AntipsychTyp == "0" & AntipsychAtyp == "0" ~ "none",
      TRUE ~ NA_character_
    ),
    neuroleptic = factor(neuroleptic, levels = c("none","typ","atyp","both")))

donor_neuroleptic <- Ruz_MtSinai@meta.data %>%
  distinct(Donor = Donor, neuroleptic)  

  #Maybe keep ruzicka as one df for this because of sample size issuses 