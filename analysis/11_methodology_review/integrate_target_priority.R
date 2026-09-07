options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

base <- "D:/PDAC_P1"
out <- file.path(base, "analysis/11_methodology_review")
priority <- fread(file.path(base, "analysis/08_program_validation/endothelial_target_prioritization.tsv"), sep = "\t", na.strings = c("", "NA"))
coverage <- fread(file.path(base, "analysis/10_communication_pseudotime/Target_screen/endothelial_target_coverage_summary.tsv"), sep = "\t", na.strings = c("", "NA"))
dep <- fread(file.path(out, "DepMap_PDAC_candidate_dependency_summary.tsv"), sep = "\t", na.strings = c("", "NA"))
hpa <- fread(file.path(base, "analysis/10_communication_pseudotime/Target_screen/HPA_target_summary.tsv"), sep = "\t", na.strings = c("", "NA"))

z <- merge(priority, coverage[, .(gene_id = gene, mean_detection, mean_patients_expressed)], by = "gene_id", all = TRUE)
z <- merge(z, dep[, .(gene_id, depmap_n_models = n_models, depmap_median_effect = median_gene_effect,
                      depmap_frac_leq_neg05 = frac_dependency_leq_neg05)], by = "gene_id", all = TRUE)
if ("gene_id" %in% names(hpa)) z <- merge(z, hpa[, .(gene_id = gene, hpa_evidence, protein_class, cancer_specificity)], by = "gene_id", all = TRUE)

rank01 <- function(x, decreasing = TRUE) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  if (decreasing) r <- max(r, na.rm = TRUE) - r + 1
  r / max(r, na.rm = TRUE)
}
z <- z[!is.na(mean_detection) & is.finite(evidence_score)]
z[, expr_replication_score := rank01(evidence_score)]
z[, coverage_score := rank01(mean_detection)]
z[, depmap_dependency_score := rank01(-depmap_median_effect)]
z[, integrated_priority_score := 0.45 * expr_replication_score + 0.30 * coverage_score +
      0.25 * fifelse(is.na(depmap_median_effect), 0, depmap_dependency_score)]
setorder(z, -integrated_priority_score)
z[, integrated_rank := seq_len(.N)]
write.table(z, file.path(out, "PDAC_endothelial_integrated_target_priority.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

report <- c(
  "# Integrated endothelial target priority",
  "",
  "The score is a transparent triage aid, not a statistical significance test. It combines cross-cohort expression evidence (45%), single-cell patient/cell coverage (30%), and DepMap PDAC dependency (25%).",
  "",
  paste0("Top candidates: ", paste(head(z$gene_id, 5), collapse = ", "), "."),
  "Candidates with missing DepMap values are retained but cannot be claimed as functionally supported.",
  "The next gate is cis-eQTL/pQTL availability and colocalization feasibility before any MR claim.",
  "",
  "Output: PDAC_endothelial_integrated_target_priority.tsv"
)
writeLines(report, file.path(out, "PDAC_endothelial_integrated_target_priority_report.md"))
