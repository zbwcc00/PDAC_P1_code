library(Seurat)
library(Matrix)

datasets <- c(
  GSE154778 = "D:/第二篇大论文/analysis/03_scrna_annotation/seurat_annotated/GSE154778__primary_tumor.rds",
  GSE155698 = "D:/第二篇大论文/analysis/03_scrna_annotation/seurat_annotated/GSE155698__primary_tumor.rds",
  GSE212966 = "D:/第二篇大论文/analysis/03_scrna_annotation/seurat_annotated/GSE212966__primary_tumor.rds"
)
out <- "C:/Users/HUAWEI/Documents/ChatGPT/Try/vm_exploratory_output"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# VM_core contains genes recurrently used in mechanistic VM studies, while the
# endothelial_anchor set is kept separate to test whether any VM-like signal is
# merely a proxy for endothelial abundance. This is an exploratory score, not a
# validated PDAC VM classifier.
vm_core <- c("EPHA2", "CDH5", "MCAM", "LAMC2", "LAMA4", "MMP2", "MMP9",
             "ITGB1", "RHOA", "ROCK1", "ROCK2")
endothelial_anchor <- c("PECAM1", "VWF", "EMCN", "ENG", "KDR", "FLT1")

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

summaries <- list()
cell_files <- character()

for (dataset in names(datasets)) {
  object <- readRDS(datasets[[dataset]])
  metadata <- object[[]]
  assay_name <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
  expr <- GetAssayData(object, assay = assay_name, layer = "data")
  genes <- rownames(expr)
  vm_use <- intersect(vm_core, genes)
  anchor_use <- intersect(endothelial_anchor, genes)
  if (length(vm_use) < 6 || length(anchor_use) < 4) {
    stop(dataset, ": insufficient genes for exploratory score")
  }
  vm_matrix <- as.matrix(expr[vm_use, , drop = FALSE])
  anchor_matrix <- as.matrix(expr[anchor_use, , drop = FALSE])
  # Standardize each gene across cells, then average genes for each cell.
  vm_score <- colMeans(t(apply(vm_matrix, 1, zscore)))
  anchor_score <- colMeans(t(apply(anchor_matrix, 1, zscore)))
  epas1_score <- if ("EPAS1" %in% genes) as.numeric(expr["EPAS1", ]) else rep(NA_real_, ncol(expr))
  metadata$VM_core_score <- as.numeric(vm_score)
  metadata$endothelial_anchor_score <- as.numeric(anchor_score)
  metadata$EPAS1_expression <- epas1_score
  metadata$VM_core_score_resid_endothelial <- NA_real_

  # Residualize the VM-like score on endothelial anchor within each lineage,
  # avoiding a cross-lineage mixture interpretation.
  for (lin in unique(metadata$major_lineage)) {
    idx <- which(metadata$major_lineage == lin)
    if (length(idx) >= 20 && is.finite(sd(metadata$endothelial_anchor_score[idx], na.rm = TRUE)) && sd(metadata$endothelial_anchor_score[idx], na.rm = TRUE) > 0) {
      fit <- lm(VM_core_score ~ endothelial_anchor_score, data = metadata[idx, , drop = FALSE])
      metadata$VM_core_score_resid_endothelial[idx] <- resid(fit)
    }
  }

  cell_path <- file.path(out, paste0(dataset, "_cell_scores.tsv"))
  write.table(cbind(cell_id = rownames(metadata), metadata[, c("dataset", "patient_id", "major_lineage", "VM_core_score", "endothelial_anchor_score", "VM_core_score_resid_endothelial", "EPAS1_expression")]),
              cell_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cell_files <- c(cell_files, cell_path)

  keep <- metadata$major_lineage %in% c("malignant_epithelial", "fibroblast_CAF", "endothelial", "epithelial")
  sub <- metadata[keep, c("patient_id", "major_lineage", "VM_core_score", "endothelial_anchor_score", "VM_core_score_resid_endothelial", "EPAS1_expression"), drop = FALSE]
  patient_lineage <- aggregate(sub[, c("VM_core_score", "endothelial_anchor_score", "VM_core_score_resid_endothelial", "EPAS1_expression")],
                               by = sub[, c("patient_id", "major_lineage")], FUN = function(x) mean(x, na.rm = TRUE))
  patient_lineage$dataset <- dataset
  summaries[[dataset]] <- patient_lineage
}

summary_df <- do.call(rbind, summaries)
write.table(summary_df, file.path(out, "vm_patient_lineage_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

stats <- list()
for (dataset in unique(summary_df$dataset)) {
  d <- summary_df[summary_df$dataset == dataset & summary_df$major_lineage %in% c("malignant_epithelial", "fibroblast_CAF", "endothelial"), , drop = FALSE]
  for (lin in unique(d$major_lineage)) {
    x <- d[d$major_lineage == lin, , drop = FALSE]
    stats[[length(stats) + 1]] <- data.frame(
      dataset = dataset,
      lineage = lin,
      n_patients = nrow(x),
      median_VM_core = median(x$VM_core_score, na.rm = TRUE),
      median_VM_residual = median(x$VM_core_score_resid_endothelial, na.rm = TRUE),
      median_endothelial_anchor = median(x$endothelial_anchor_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}
write.table(do.call(rbind, stats), file.path(out, "vm_patient_lineage_summary_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "Exploratory lineage-resolved VM score audit",
  "VM_core genes:", paste(vm_use, collapse = ", "),
  "Endothelial anchor genes:", paste(anchor_use, collapse = ", "),
  "Interpretation: a VM-like expression score in malignant epithelial or CAF cells is not equivalent to histologically verified VM.",
  "Required orthogonal validation: PAS+/CD31- channels containing erythrocytes and malignant/CAF lineage markers, ideally with spatial protein imaging."
), file.path(out, "vm_score_audit_notes.txt"))
