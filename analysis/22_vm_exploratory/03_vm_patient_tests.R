out <- "C:/Users/HUAWEI/Documents/ChatGPT/Try/vm_exploratory_output"
summary_df <- read.table(file.path(out, "vm_patient_lineage_summary.tsv"), sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

lineages <- c("endothelial", "fibroblast_CAF", "malignant_epithelial")
tests <- list()
for (dataset in unique(summary_df$dataset)) {
  d <- summary_df[summary_df$dataset == dataset & summary_df$major_lineage %in% lineages, , drop = FALSE]
  for (pair in list(c("endothelial", "fibroblast_CAF"), c("endothelial", "malignant_epithelial"), c("fibroblast_CAF", "malignant_epithelial"))) {
    a <- d[d$major_lineage == pair[1], "VM_core_score"]
    b <- d[d$major_lineage == pair[2], "VM_core_score"]
    tests[[length(tests) + 1]] <- data.frame(
      dataset = dataset,
      comparison = paste(pair, collapse = " vs "),
      n_a = length(a), n_b = length(b),
      median_a = median(a, na.rm = TRUE), median_b = median(b, na.rm = TRUE),
      wilcox_p = wilcox.test(a, b, exact = FALSE)$p.value,
      stringsAsFactors = FALSE
    )
  }
}
tests_df <- do.call(rbind, tests)
tests_df$BH_q <- p.adjust(tests_df$wilcox_p, method = "BH")
write.table(tests_df, file.path(out, "vm_patient_lineage_pairwise_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Cross-compartment coupling: compare patient-level malignant/CAF residual VM
# scores with endothelial EPAS1 and endothelial-anchor scores. This is a
# hypothesis-generating association, not a test of VM formation.
coupling <- list()
for (dataset in unique(summary_df$dataset)) {
  d <- summary_df[summary_df$dataset == dataset, , drop = FALSE]
  endo <- d[d$major_lineage == "endothelial", c("patient_id", "VM_core_score", "VM_core_score_resid_endothelial", "endothelial_anchor_score")]
  for (lin in c("fibroblast_CAF", "malignant_epithelial")) {
    other <- d[d$major_lineage == lin, c("patient_id", "VM_core_score", "VM_core_score_resid_endothelial")]
    merged <- merge(endo, other, by = "patient_id", suffixes = c("_endothelial", "_other"))
    if (nrow(merged) >= 4) {
      r1 <- suppressWarnings(cor.test(merged$VM_core_score_resid_endothelial_other, merged$VM_core_score_endothelial, method = "spearman", exact = FALSE))
      r2 <- suppressWarnings(cor.test(merged$VM_core_score_resid_endothelial_other, merged$endothelial_anchor_score, method = "spearman", exact = FALSE))
      epas1 <- d[d$major_lineage == "endothelial", c("patient_id", "EPAS1_expression")]
      merged_epas1 <- merge(merged, epas1, by = "patient_id")
      epas1_ok <- is.finite(merged_epas1$EPAS1_expression) & is.finite(merged_epas1$VM_core_score_resid_endothelial_other)
      if (sum(epas1_ok) >= 4 && length(unique(merged_epas1$EPAS1_expression[epas1_ok])) > 1 && length(unique(merged_epas1$VM_core_score_resid_endothelial_other[epas1_ok])) > 1) {
        r3 <- suppressWarnings(cor.test(merged_epas1$VM_core_score_resid_endothelial_other[epas1_ok], merged_epas1$EPAS1_expression[epas1_ok], method = "spearman", exact = FALSE))
        rho_epas1 <- unname(r3$estimate)
        p_epas1 <- r3$p.value
      } else {
        rho_epas1 <- NA_real_
        p_epas1 <- NA_real_
      }
      coupling[[length(coupling) + 1]] <- data.frame(
        dataset = dataset, other_lineage = lin, n_patients = nrow(merged),
        rho_other_residual_vs_endothelial_VM = unname(r1$estimate),
        p_other_residual_vs_endothelial_VM = r1$p.value,
        rho_other_residual_vs_endothelial_anchor = unname(r2$estimate),
        p_other_residual_vs_endothelial_anchor = r2$p.value,
        rho_other_residual_vs_endothelial_EPAS1 = rho_epas1,
        p_other_residual_vs_endothelial_EPAS1 = p_epas1,
        stringsAsFactors = FALSE
      )
    }
  }
}
coupling_df <- do.call(rbind, coupling)
coupling_df$BH_q_EPAS1 <- p.adjust(coupling_df$p_other_residual_vs_endothelial_EPAS1, method = "BH")
write.table(coupling_df, file.path(out, "vm_cross_compartment_coupling_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "Interpretive rule:",
  "1) A higher VM-like expression score in endothelial cells is expected because the score includes shared vascular genes.",
  "2) A malignant/CAF score is not histological VM without PAS+/CD31- channels, erythrocytes, and lineage markers.",
  "3) Cross-compartment correlations are exploratory and do not establish endothelial-to-tumor conversion or perfusable VM."
), file.path(out, "vm_patient_test_interpretation.txt"))
