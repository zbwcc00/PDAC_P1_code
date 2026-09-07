options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

root <- "D:/PDAC_P1"
outdir <- file.path(root, "analysis/17_P0_reinforcement")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# The two discovery/validation cohorts share the only defensible phenotype contrast:
# primary PDAC versus adjacent normal, within endothelial pseudobulk samples.
de_root <- file.path(root, "analysis/06_pseudobulk_DE")
cohorts <- c("GSE155698", "GSE212966")

read_cohort <- function(dataset) {
  edge <- fread(file.path(de_root, "edgeR", dataset, "endothelial_DE.tsv"))[
    , .(gene_id, edge_logFC = logFC, edge_FDR = FDR)]
  deseq <- fread(file.path(de_root, "DESeq2", dataset, "endothelial_DE.tsv"))[
    , .(gene_id, deseq_logFC = log2FoldChange, deseq_FDR = padj)]
  merge(edge, deseq, by = "gene_id")
}

dat <- setNames(lapply(cohorts, read_cohort), cohorts)

# Pre-specification: a discovery score only rewards within-cohort agreement between
# edgeR and DESeq2. No spatial information is used for ranking.
score_discovery <- function(x) {
  x <- copy(x)
  x[, method_direction_agreement := is.finite(edge_logFC) & is.finite(deseq_logFC) &
      sign(edge_logFC) == sign(deseq_logFC) & edge_logFC != 0 & deseq_logFC != 0]
  x[, n_method_FDR05 := as.integer(edge_FDR < 0.05) + as.integer(deseq_FDR < 0.05)]
  x[, median_logFC := apply(as.matrix(.SD), 1, median, na.rm = TRUE), .SDcols = c("edge_logFC", "deseq_logFC")]
  x[, mean_abs_logFC := rowMeans(abs(as.matrix(.SD)), na.rm = TRUE), .SDcols = c("edge_logFC", "deseq_logFC")]
  setorder(x, -method_direction_agreement, -n_method_FDR05, -mean_abs_logFC, edge_FDR, deseq_FDR)
  x[, discovery_rank := seq_len(.N)]
  x
}

loco <- rbindlist(lapply(cohorts, function(held_out) {
  discovery <- setdiff(cohorts, held_out)
  discovery_tab <- score_discovery(dat[[discovery]])
  validation <- copy(dat[[held_out]])
  setnames(validation, setdiff(names(validation), "gene_id"), paste0("heldout_", setdiff(names(validation), "gene_id")))
  result <- merge(discovery_tab, validation, by = "gene_id", all.x = TRUE)
  result[, c("discovery_cohort", "heldout_cohort") := list(discovery, held_out)]
  result[, heldout_direction_match := is.finite(heldout_edge_logFC) & is.finite(heldout_deseq_logFC) &
      sign(heldout_edge_logFC) == sign(median_logFC) & sign(heldout_deseq_logFC) == sign(median_logFC)]
  result[, heldout_any_FDR05 := (heldout_edge_FDR < 0.05 | heldout_deseq_FDR < 0.05)]
  result[, heldout_both_FDR05 := (heldout_edge_FDR < 0.05 & heldout_deseq_FDR < 0.05)]
  result
}), fill = TRUE)

fwrite(loco, file.path(outdir, "loco_endothelial_all_gene_ranks.tsv"), sep = "\t", na = "NA")
epas1_loco <- loco[gene_id == "EPAS1"]
fwrite(epas1_loco, file.path(outdir, "loco_EPAS1_summary.tsv"), sep = "\t", na = "NA")

# GSE154778 does not contain an endothelial pseudobulk contrast and contrasts
# primary tumors against metastases; it is therefore recorded as non-eligible,
# rather than treated as an invalid third LOCO fold.
loco_scope <- data.table(
  cohort = c("GSE155698", "GSE212966", "GSE154778"),
  endothelial_DE_available = c(TRUE, TRUE, FALSE),
  contrast = c("primary_PDAC_vs_adjacent_normal", "primary_PDAC_vs_adjacent_normal", "primary_vs_metastasis"),
  loco_eligible = c(TRUE, TRUE, FALSE),
  rationale = c(
    "shared phenotype contrast", "shared phenotype contrast",
    "No endothelial pseudobulk contrast; different clinical comparison"
  )
)
fwrite(loco_scope, file.path(outdir, "loco_scope_and_eligibility.tsv"), sep = "\t")

# Patient-level spatial sensitivity. Spots estimate each patient's correlation;
# the meta-analysis unit remains patient/sample, not individual spots.
spatial_root <- file.path(root, "analysis/09_spatial_validation")
spatial_sets <- c("GSE282302", "GSE297144")

cor_spearman <- function(a, b) {
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < 10L || length(unique(a[keep])) < 2L || length(unique(b[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor(a[keep], b[keep], method = "spearman"))
}

spatial_patient <- rbindlist(lapply(spatial_sets, function(dataset) {
  x <- fread(file.path(spatial_root, dataset, paste0(dataset, "_spot_anchor_scores.tsv")))
  x <- x[in_tissue == 1]
  x[, patient_id := as.character(patient)]
  x[, batch_id := as.character(sample)]
  rbindlist(lapply(split(x, x$patient_id), function(d) {
    covars <- data.frame(log_total_counts = log1p(d$total_counts), log_n_genes = log1p(d$n_genes))
    epas1_resid <- residuals(lm(EPAS1_logCPM ~ log_total_counts + log_n_genes, data = cbind(d, covars)))
    vascular_resid <- residuals(lm(vascular_marker_score ~ log_total_counts + log_n_genes, data = cbind(d, covars)))
    rho_adjusted <- cor_spearman(epas1_resid, vascular_resid)
    data.table(
      dataset = dataset,
      patient = d$patient_id[1],
      n_spots = nrow(d),
      n_sections = uniqueN(d$batch_id),
      rho_unadjusted = cor_spearman(d$EPAS1_logCPM, d$vascular_marker_score),
      rho_adjusted_depth_detection = rho_adjusted,
      fisher_z_adjusted = atanh(pmin(pmax(rho_adjusted, -0.999999), 0.999999)),
      var_fisher_z = ifelse(nrow(d) > 3, 1 / (nrow(d) - 3), NA_real_)
    )
  }))
}), fill = TRUE)
fwrite(spatial_patient, file.path(outdir, "spatial_patient_level_adjusted_correlations.tsv"), sep = "\t", na = "NA")

random_effects <- function(d) {
  d <- d[is.finite(fisher_z_adjusted) & is.finite(var_fisher_z) & var_fisher_z > 0]
  if (nrow(d) < 2) return(data.table(k = nrow(d), pooled_rho_RE = NA_real_, ci_low_rho_RE = NA_real_, ci_high_rho_RE = NA_real_, tau2_DL = NA_real_, I2_percent = NA_real_))
  w_fixed <- 1 / d$var_fisher_z
  z_fixed <- sum(w_fixed * d$fisher_z_adjusted) / sum(w_fixed)
  q <- sum(w_fixed * (d$fisher_z_adjusted - z_fixed)^2)
  c_value <- sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)
  tau2 <- max(0, (q - (nrow(d) - 1)) / c_value)
  w_re <- 1 / (d$var_fisher_z + tau2)
  z_re <- sum(w_re * d$fisher_z_adjusted) / sum(w_re)
  se_re <- sqrt(1 / sum(w_re))
  data.table(k = nrow(d), pooled_rho_RE = tanh(z_re), ci_low_rho_RE = tanh(z_re - 1.96 * se_re), ci_high_rho_RE = tanh(z_re + 1.96 * se_re), tau2_DL = tau2, I2_percent = ifelse(q > 0, max(0, (q - (nrow(d) - 1)) / q) * 100, 0))
}

spatial_meta <- rbindlist(c(
  list(cbind(scope = "all_independent_patients_or_samples", random_effects(spatial_patient))),
  lapply(spatial_sets, function(dset) cbind(scope = dset, random_effects(spatial_patient[dataset == dset])))
), fill = TRUE)
fwrite(spatial_meta, file.path(outdir, "spatial_adjusted_correlation_random_effects_meta.tsv"), sep = "\t", na = "NA")

spatial_scope <- spatial_patient[, .(n_patient_units = .N, n_multi_section_units = sum(n_sections > 1), max_sections_per_unit = max(n_sections)), by = dataset]
spatial_scope[, independence_note := ifelse(dataset == "GSE297144", "One ROI/sample per patient label; donor identifiers unavailable in processed GEO metadata, so analysis unit is sample.", "Multiple ROIs aggregated within each patient before meta-analysis.")]
fwrite(spatial_scope, file.path(outdir, "spatial_independence_audit.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(outdir, "01_sessionInfo.txt"))
