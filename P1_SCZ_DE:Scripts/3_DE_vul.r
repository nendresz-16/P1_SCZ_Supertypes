# run in conda activate de_env
setwd("P1_SCZ_DE")

library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(edgeR)
library(limma)
library(EnhancedVolcano)
library(patchwork)


meta_all <- read.csv("Files/pseudobulk_metadata_SCZ.csv", check.names = FALSE)

# 1 row per Donor × Vulnerability
meta_all2 <- meta_all %>%
  group_by(Donor, Vulnerability) %>%
  summarise(
    Diagnosis = first(Diagnosis),
    Sex       = first(Sex),
    Age       = first(Age),
    PMI       = first(PMI),
    num_cells = sum(num_cells),
    .groups   = "drop"
  )

# sanity: (Donor, Vulnerability) must be unique
stopifnot(nrow(meta_all2) == nrow(distinct(meta_all2, Donor, Vulnerability)))

# List cohorts explicitly (safer than a broad pattern)
pb_files <- c(
  "Files/Pseudobulk_OFC.rds",
  "Files/Pseudobulk_Bat.rds",
  "Files/Pseudobulk_McLean.rds",
  "Files/Pseudobulk_MtSinai.rds"
)
pb_files <- pb_files[file.exists(pb_files)]
if (length(pb_files) == 0) stop("No pseudobulk files found in Files/.")

run_within_diagnosis <- function(pb_path, cohort_label,
                                 min_prop = 0.80, min_per_group = 4,
                                 fdr_cut = 0.05) {
  message("\n=== ", cohort_label, " ===")
  counts <- readRDS(pb_path)
  colnames(counts) <- gsub("-", "_", colnames(counts))  # normalize dashes

  # ---- build sample metadata (robust split) ----
  smeta <- tibble(SampleID = colnames(counts)) %>%
    # split on the FIRST underscore only; keeps donors that may contain underscores
    mutate(
      Vulnerability = sub("^(.*?)_.*$", "\\1", SampleID),
      Donor         = sub("^.*?_(.*)$", "\\1", SampleID)
    )

  # If splitting failed (no underscore), try a fallback (pipe)
  if (any(smeta$Vulnerability == smeta$SampleID | smeta$Donor == smeta$SampleID)) {
    smeta <- tibble(SampleID = colnames(counts)) %>%
      mutate(
        Vulnerability = sub("^(.*?)\\|.*$", "\\1", SampleID),
        Donor         = sub("^.*?\\|(.*)$", "\\1", SampleID)
      )
  }

  # Validate split
  bad <- which(is.na(smeta$Vulnerability) | is.na(smeta$Donor) |
                 smeta$Vulnerability == "" | smeta$Donor == "")
  if (length(bad)) {
    warning(cohort_label, ": could not parse some SampleIDs. Examples: ",
            paste(smeta$SampleID[bad][1:min(5, length(bad))], collapse = ", "))
    return(invisible(NULL))
  }

  smeta <- smeta %>%
    left_join(meta_all2, by = c("Donor","Vulnerability"))

  # If join failed for some columns, report and skip
  na_rows <- which(is.na(smeta$Diagnosis) | is.na(smeta$num_cells))
  if (length(na_rows)) {
    warning(cohort_label, ": unmatched (Donor,Vulnerability) for ",
            length(na_rows), " samples. Examples: ",
            paste(paste0(smeta$Vulnerability[na_rows][1:min(5,length(na_rows))], "_",
                         smeta$Donor[na_rows][1:min(5,length(na_rows))]), collapse = ", "))
    return(invisible(NULL))
  }

  smeta <- smeta %>%
    mutate(
      Vulnerability = factor(Vulnerability, levels = c("Unaffected","Affected")),
      Diagnosis     = factor(Diagnosis,   levels = c("Control","Schizophrenia")),
      log_10_cells_per_donor = log10(num_cells)
    )

  # Align counts ↔ smeta
  counts <- counts[, smeta$SampleID, drop = FALSE]
  if (!identical(colnames(counts), smeta$SampleID)) {
    warning(cohort_label, ": counts colnames and smeta$SampleID not identical; skipping.")
    return(invisible(NULL))
  }

  message("Sample table (Dx × Vuln):")
  print(with(smeta, table(Diagnosis, Vulnerability)))

  # ---- DE: Affected vs Unaffected within each Dx ----
  DE_out <- list(); plots <- list()
  for (dx in c("Control","Schizophrenia")) {
    smeta_dx <- smeta %>% filter(Diagnosis == dx)
    n_aff  <- sum(smeta_dx$Vulnerability == "Affected",   na.rm = TRUE)
    n_unaf <- sum(smeta_dx$Vulnerability == "Unaffected", na.rm = TRUE)
    message("  ", dx, ": nUnaff=", n_unaf, " nAff=", n_aff)

    if (n_aff < min_per_group || n_unaf < min_per_group) {
      warning("Skipping ", cohort_label, " within ", dx, " (not enough samples).")
      next
    }

    counts_dx <- counts[, smeta_dx$SampleID, drop = FALSE]

    dge <- edgeR::DGEList(counts = counts_dx, genes = rownames(counts_dx))
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

    min_samples <- ceiling(ncol(dge) * min_prop)
    keep_genes  <- rowSums(dge$counts >= 1) >= min_samples
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

    design <- model.matrix(~ Age + PMI + Sex + log_10_cells_per_donor + Vulnerability, data = smeta_dx)
    vm  <- limma::voom(dge, design, plot = FALSE)
    fit <- limma::eBayes(limma::lmFit(vm, design))

    DE <- limma::topTable(fit, coef = "VulnerabilityAffected", n = Inf,
                          adjust.method = "BH", sort = "none")
    DE_out[[dx]] <- DE

    sig_n <- sum(DE$adj.P.Val < fdr_cut, na.rm = TRUE)
    plots[[dx]] <- EnhancedVolcano(
      DE,
      lab = rownames(DE),
      x = "logFC",
      y = "adj.P.Val",
      title = paste0(cohort_label, " | ", dx, " | Aff vs Unaff\nFDR<", fdr_cut, ": ", sig_n,
                     " | nUnaff=", n_unaf, " nAff=", n_aff),
      pCutoff = fdr_cut,
      FCcutoff = 1.0,
      pointSize = 2.0,
      labSize = 3.0
    )
  }

  # Save combined results and grid (even if only one side exists)
  saveRDS(DE_out, file.path("Files", paste0("DE_", cohort_label, "_Vuln_withinDx.rds")))
  if (length(plots)) {
    combined <- if (all(c("Control","Schizophrenia") %in% names(plots))) {
      plots$Control + plots$Schizophrenia + patchwork::plot_layout(ncol = 2)
    } else {
      patchwork::wrap_plots(plots)
    }
    ggsave(file.path("Figures", paste0("Volcano_", cohort_label, "_Vuln_withinDx_grid.png")),
           combined, width = 14, height = 6, dpi = 300)
  }
  TRUE
}

# ---- Run over all cohorts with error isolation ----
for (fp in pb_files) {
  cohort <- sub("^Pseudobulk_(.*)\\.rds$", "\\1", basename(fp))
  tryCatch(
    run_within_diagnosis(fp, cohort_label = cohort, fdr_cut = 0.05),
    error = function(e) {
      message("ERROR in ", cohort, ": ", conditionMessage(e))
    }
  )
}

# See what was produced
print(list.files("Files",   pattern = "^DE_.*_Vuln_withinDx\\.rds$", full.names = TRUE))
print(list.files("Figures", pattern = "^Volcano_.*_Vuln_withinDx_grid\\.png$", full.names = TRUE))


#for batiuk only 3 in each group 


meta_all <- read.csv("Files/pseudobulk_metadata_SCZ.csv", check.names = FALSE)
meta_all2 <- meta_all %>%
  group_by(Donor, Vulnerability) %>%
  summarise(
    Diagnosis = first(Diagnosis),
    Sex = first(Sex), Age = first(Age), PMI = first(PMI),
    num_cells = sum(num_cells), .groups = "drop"
  )

pb_path <- "Files/Pseudobulk_Bat.rds"
counts  <- readRDS(pb_path)
colnames(counts) <- gsub("-", "_", colnames(counts))

smeta <- tibble(SampleID = colnames(counts)) %>%
  mutate(
    Vulnerability = sub("^(.*?)_.*$", "\\1", SampleID),
    Donor         = sub("^.*?_(.*)$", "\\1", SampleID)
  ) %>%
  left_join(meta_all2, by = c("Donor","Vulnerability")) %>%
  mutate(
    Vulnerability = factor(Vulnerability, levels = c("Unaffected","Affected")),
    Diagnosis     = factor(Diagnosis,   levels = c("Control","Schizophrenia")),
    log_10_cells_per_donor = log10(num_cells)
  )

counts <- counts[, smeta$SampleID, drop = FALSE]

DE_out <- list(); plots <- list()
min_prop <- 0.80; min_per_group <- 3; fdr_cut <- 0.05

for (dx in c("Control","Schizophrenia")) {
  smeta_dx <- smeta %>% filter(Diagnosis == dx)
  n_aff  <- sum(smeta_dx$Vulnerability == "Affected")
  n_unaf <- sum(smeta_dx$Vulnerability == "Unaffected")
  if (n_aff < min_per_group || n_unaf < min_per_group) next

  counts_dx <- counts[, smeta_dx$SampleID, drop = FALSE]
  dge <- edgeR::DGEList(counts = counts_dx, genes = rownames(counts_dx))
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  keep <- rowSums(dge$counts >= 1) >= ceiling(ncol(dge) * min_prop)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  design <- model.matrix(~ Age + PMI + Sex + log_10_cells_per_donor + Vulnerability, data = smeta_dx)
  vm <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::eBayes(limma::lmFit(vm, design))
  DE <- limma::topTable(fit, coef = "VulnerabilityAffected", n = Inf, adjust.method = "BH", sort = "none")
  DE_out[[dx]] <- DE

  sig_n <- sum(DE$adj.P.Val < fdr_cut, na.rm = TRUE)
  plots[[dx]] <- EnhancedVolcano(
    DE, lab = rownames(DE), x = "logFC", y = "adj.P.Val",
    title = paste0("Bat | ", dx, " | Aff vs Unaff\nFDR<", fdr_cut, ": ", sig_n,
                   " | nUnaff=", n_unaf, " nAff=", n_aff),
    pCutoff = fdr_cut, FCcutoff = 1.0, pointSize = 2.0, labSize = 3.0
  )
}

saveRDS(DE_out, "Files/DE_Bat_Vuln_withinDx.rds")
if (length(plots)) {
  combined <- if (all(c("Control","Schizophrenia") %in% names(plots))) {
    plots$Control + plots$Schizophrenia + patchwork::plot_layout(ncol = 2)
  } else wrap_plots(plots)
  ggsave("Figures/Volcano_Bat_Vuln_withinDx_grid.png", combined, width = 14, height = 6, dpi = 300)
}