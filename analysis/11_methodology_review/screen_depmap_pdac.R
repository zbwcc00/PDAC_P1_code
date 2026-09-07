options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

base <- "D:/PDAC_P1"
raw <- file.path(base, "data/03_functional/DepMap24Q4/raw")
out <- file.path(base, "analysis/11_methodology_review")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

model_file <- file.path(raw, "Model.csv")
crispr_file <- file.path(raw, "CRISPRGeneEffect.csv")
priority_file <- file.path(base, "analysis/08_program_validation/endothelial_target_prioritization.tsv")

model <- fread(model_file, na.strings = c("", "NA"), showProgress = FALSE)
priority <- fread(priority_file, sep = "\t", showProgress = FALSE)
genes <- unique(priority$gene_id)
genes <- genes[!is.na(genes) & nzchar(genes)]

pdac_mask <- grepl("^Pancreas$", model$OncotreeLineage, ignore.case = TRUE) &
  grepl("Pancreatic Adenocarcinoma|Adenosquamous Carcinoma of the Pancreas",
        model$OncotreePrimaryDisease, ignore.case = TRUE)
pdac_models <- model[pdac_mask]

header <- names(fread(crispr_file, nrows = 0, showProgress = FALSE))
gene_cols <- setNames(vapply(genes, function(g) {
  hit <- header[grepl(paste0("^", g, " \\("), header)]
  if (length(hit)) hit[1] else NA_character_
}, character(1)), genes)
gene_cols <- gene_cols[!is.na(gene_cols)]

keep_ids <- intersect(pdac_models$ModelID, sub("^X", "", header[1]))
crispr <- fread(crispr_file, select = c(header[1], unname(gene_cols)), na.strings = c("", "NA"), showProgress = FALSE)
setnames(crispr, 1, "ModelID")
crispr$ModelID <- sub("^X", "", crispr$ModelID)
pdac_crispr <- crispr[ModelID %in% pdac_models$ModelID]

long <- rbindlist(lapply(names(gene_cols), function(g) {
  x <- suppressWarnings(as.numeric(pdac_crispr[[gene_cols[[g]]]]))
  data.table(gene_id = g, n_models = sum(!is.na(x)), median_gene_effect = median(x, na.rm = TRUE),
             mean_gene_effect = mean(x, na.rm = TRUE), sd_gene_effect = sd(x, na.rm = TRUE),
             frac_strong_dependency_leq_neg1 = mean(x <= -1, na.rm = TRUE),
             frac_dependency_leq_neg05 = mean(x <= -0.5, na.rm = TRUE))
}), fill = TRUE)
setorder(long, median_gene_effect)

annot <- pdac_models[, .(ModelID, CellLineName, StrippedCellLineName, OncotreeLineage,
                         OncotreePrimaryDisease, OncotreeSubtype, OncotreeCode,
                         PrimaryOrMetastasis, SourceType)]
write.table(long, file.path(out, "DepMap_PDAC_candidate_dependency_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annot, file.path(out, "DepMap_PDAC_model_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

report <- c(
  "# DepMap PDAC dependency screen",
  "",
  paste0("- DepMap release: 24Q4 public; model rows: ", nrow(model), "; PDAC-like models selected: ", nrow(pdac_models), "."),
  paste0("- Candidate genes requested: ", length(genes), "; genes available in CRISPRGeneEffect: ", length(gene_cols), "."),
  "- Gene effect scores are interpreted as dependency when more negative; summaries are descriptive and require lineage-aware controls.",
  "- Strong dependency fraction thresholds (≤−1 and ≤−0.5) are reported for prioritization, not used as universal biological cutoffs.",
  "",
  "Outputs: DepMap_PDAC_candidate_dependency_summary.tsv and DepMap_PDAC_model_annotation.tsv"
)
writeLines(report, file.path(out, "DepMap_PDAC_dependency_report.md"))
