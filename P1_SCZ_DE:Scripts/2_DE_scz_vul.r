# run on conda activate de_env
setwd("P1_SCZ_DE")
library(dplyr)
library(tidyr)
library(edgeR)
library(limma)
library(tibble)
library(dplyr)
library(tidyr)
library(tibble)
library(Matrix)
library(edgeR)
library(limma)
library(EnhancedVolcano)
library(patchwork)


meta_all <- read.csv("Files/pseudobulk_metadata_SCZ.csv", check.names = FALSE)

# one row per Donor×Vulnerability
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

# list the four pseudobulk files
pb_files <- list.files("Files", pattern = "^Pseudobulk_.*\\.rds$", full.names = TRUE)

run_cohort <- function(pb_path, cohort_label) {
  message("=== ", cohort_label, " ===")

  counts <- readRDS(pb_path)                # dgCMatrix or matrix (genes × samples)
  colnames(counts) <- gsub("-", "_", colnames(counts))  # normalize names

  # build sample metadata from "Unaffected_DONOR" / "Affected_DONOR"
  smeta <- tibble(SampleID = colnames(counts)) %>%
    separate(SampleID, into = c("Vulnerability","Donor"),
             sep = "_", extra = "merge", fill = "right", remove = FALSE) %>%
    left_join(meta_all2, by = c("Donor","Vulnerability")) %>%
    mutate(
      Vulnerability = factor(Vulnerability, levels = c("Unaffected","Affected")),
      Diagnosis     = factor(Diagnosis,   levels = c("Control","Schizophrenia")),
      log_10_cells_per_donor = log10(num_cells)
    )

  # align counts to metadata rows
  counts <- counts[, smeta$SampleID, drop = FALSE]
  stopifnot(identical(colnames(counts), smeta$SampleID))

  # ---- DE per vulnerability (SCZ vs Control) ----
  DE_out <- list()
  MIN_PER_GRP <- 4
  MIN_PROP    <- 0.80

  for (vuln in levels(smeta$Vulnerability)) {
    smeta_sub <- smeta %>% filter(Vulnerability == vuln, !is.na(Diagnosis))
    n_SCZ <- sum(smeta_sub$Diagnosis == "Schizophrenia", na.rm = TRUE)
    n_CON <- sum(smeta_sub$Diagnosis == "Control",       na.rm = TRUE)
    if (n_SCZ < MIN_PER_GRP || n_CON < MIN_PER_GRP) {
      warning("Skipping ", cohort_label, " / ", vuln, " (SCZ=", n_SCZ, ", CON=", n_CON, ")")
      next
    }

    counts_sub <- counts[, smeta_sub$SampleID, drop = FALSE]

    dge <- DGEList(counts = counts_sub, genes = rownames(counts_sub))
    dge <- calcNormFactors(dge, method = "TMM")

    min_samples <- ceiling(ncol(dge) * MIN_PROP)
    keep_genes  <- rowSums(dge$counts >= 1) >= min_samples
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

    design <- model.matrix(~ Age + PMI + Sex + log_10_cells_per_donor + Diagnosis, data = smeta_sub)

    vm  <- voom(dge, design, plot = FALSE)
    fit <- eBayes(lmFit(vm, design))

    DE <- topTable(fit, coef = "DiagnosisSchizophrenia", n = Inf, adjust.method = "BH", sort = "none")
    DE$Vulnerability <- vuln
    DE$n_SCZ <- n_SCZ
    DE$n_CON <- n_CON
    DE_out[[vuln]] <- DE
  }

  # save DE
  saveRDS(DE_out, file.path("Files", paste0("DE_", cohort_label, ".rds")))

  # ---- Volcano grid (if any results) ----
  plots <- list()
  for (vuln in names(DE_out)) {
    DE <- DE_out[[vuln]]
    sig_n <- sum(DE$adj.P.Val < 0.05, na.rm = TRUE)
    ttl <- paste0(cohort_label, " ", vuln, " (SCZ vs CON)\nFDR<0.05: ", sig_n,
                  " | nSCZ=", unique(DE$n_SCZ), " nCON=", unique(DE$n_CON))

    plots[[vuln]] <- EnhancedVolcano(
      DE,
      lab = rownames(DE),
      x = "logFC",
      y = "adj.P.Val",
      title = ttl,
      pCutoff = 0.05,
      FCcutoff = 1.0,
      pointSize = 2.0,
      labSize = 3.0
    )
  }

  if (length(plots) == 0) return(invisible(NULL))

  combined <- if (all(c("Unaffected","Affected") %in% names(plots))) {
    plots$Unaffected + plots$Affected + plot_layout(ncol = 2)
  } else {
    wrap_plots(plots)
  }

  ggsave(
    filename = file.path("Figures", paste0("Volcano_", cohort_label, "_grid.png")),
    plot = combined, width = 14, height = 6, dpi = 300
  )

  invisible(TRUE)
}

# run for each cohort file
for (fp in pb_files) {
  cohort <- sub("^Pseudobulk_(.*)\\.rds$", "\\1", basename(fp))  # e.g., OFC, Bat, McLean, MtSinai
  run_cohort(fp, cohort_label = cohort)
}

# list outputs
list.files("Files", pattern = "^DE_.*\\.rds$", full.names = TRUE)
list.files("Figures", pattern = "^Volcano_.*_grid\\.png$", full.names = TRUE)
