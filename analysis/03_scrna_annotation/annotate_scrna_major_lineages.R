options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- pdac_paths()$project_root
in_dir <- file.path(root, "analysis/02_scrna_qc/seurat_checkpoints")
out_dir <- file.path(root, "analysis/03_scrna_annotation")
ann_dir <- file.path(out_dir, "seurat_annotated")
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)

marker_sets <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT17", "MUC1"),
  malignant_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT17", "MSLN", "CEACAM6", "MUC1", "S100A10"),
  fibroblast_CAF = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "SPARC", "FAP", "PDGFRA", "PDGFRB", "ACTA2"),
  myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "LILRB1", "LGALS3", "LYZ", "S100A8", "S100A9"),
  T_NK = c("CD3D", "CD3E", "TRBC1", "TRBC2", "NKG7", "GNLY", "KLRD1"),
  CD8 = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GZMK", "GZMB", "PRF1"),
  B_cell = c("CD79A", "MS4A1", "CD37", "CD74", "HLA-DRA", "CD79B"),
  endothelial = c("EMCN", "PECAM1", "VWF", "KDR", "RAMP2", "ESAM", "ENG"),
  mast = c("KIT", "TPSAB1", "MS4A2", "CPA3"),
  acinar = c("PRSS1", "REG1A", "REG1B", "AMY2A", "CPA1", "CTRB1")
)
module_sets <- list(
  malignant_epithelial_program = marker_sets$malignant_epithelial,
  CAF_program = marker_sets$fibroblast_CAF,
  myeloid_program = marker_sets$myeloid,
  CD8_activation = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GNLY", "GZMK", "GZMB", "PRF1"),
  CD8_exhaustion = c("PDCD1", "LAG3", "TIGIT", "TOX", "HAVCR2", "CTLA4", "TNFRSF9", "ENTPD1")
)

score_genes <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (!length(genes)) return(rep(0, ncol(expr)))
  Matrix::colMeans(expr[genes, , drop = FALSE])
}

classify <- function(scores) {
  best <- max.col(scores, ties.method = "first")
  best_score <- scores[cbind(seq_len(nrow(scores)), best)]
  sorted <- t(apply(scores, 1, sort, decreasing = TRUE))
  margin <- sorted[, 1] - sorted[, 2]
  labels <- colnames(scores)[best]
  labels[best_score < 0.25 | margin < 0.08] <- "unknown_ambiguous"
  labels
}

files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
cell_rows <- list(); sample_rows <- list(); errors <- list()
for (f in files) {
  key <- sub("\\.rds$", "", basename(f))
  message("Annotating ", key)
  obj <- tryCatch(readRDS(f), error = function(e) { errors <<- append(errors, list(data.frame(file = f, error = conditionMessage(e)))); NULL })
  if (is.null(obj)) next
  if ("RNA" %in% names(obj@assays)) obj <- JoinLayers(obj, assay = "RNA")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")
  lineage_scores <- sapply(marker_sets, function(g) score_genes(expr, g))
  lineage_scores <- as.matrix(lineage_scores); colnames(lineage_scores) <- names(marker_sets)
  obj$major_lineage <- classify(lineage_scores)
  for (j in seq_len(ncol(lineage_scores))) obj[[paste0("score_", colnames(lineage_scores)[j])]] <- lineage_scores[, j]
  for (nm in names(module_sets)) obj[[paste0("module_", nm)]] <- score_genes(expr, module_sets[[nm]])
  obj@misc$annotation_policy <- list(method = "marker score with conservative max-score/margin rule", marker_sets = marker_sets,
                                      unknown_rule = "max score <0.25 or top-two margin <0.08", raw_counts_preserved = TRUE)
  saveRDS(obj, file.path(ann_dir, paste0(key, ".rds")), compress = "xz")
  md <- as.data.frame(obj[[]])
  md$cell <- rownames(md)
  md$sample_id <- ifelse(!is.na(md$gsm), md$gsm, sub(".*_", "", md$orig.ident))
  md$pass_singlet <- md$passes_core_QC %in% TRUE & (is.na(md$doublet_class) | md$doublet_class != "doublet")
  keep <- md$pass_singlet
  if (any(keep)) {
    cr <- as.data.table(md[keep, c("dataset", "gsm", "tissue_context", "patient_id", "major_lineage", "module_malignant_epithelial_program", "module_CAF_program", "module_myeloid_program", "module_CD8_activation", "module_CD8_exhaustion"), drop = FALSE])
    cell_rows <<- append(cell_rows, list(cr))
    props <- cr[, .(cells = .N,
                    malignant_epithelial_score = mean(module_malignant_epithelial_program),
                    CAF_score = mean(module_CAF_program), myeloid_score = mean(module_myeloid_program),
                    CD8_activation_score = mean(module_CD8_activation), CD8_exhaustion_score = mean(module_CD8_exhaustion)),
                 by = .(dataset, gsm, tissue_context, patient_id, major_lineage)]
    sample_rows <<- append(sample_rows, list(props))
  }
  rm(obj, expr); gc(verbose = FALSE)
}

cells <- rbindlist(cell_rows, fill = TRUE)
fwrite(cells, file.path(out_dir, "PDAC_scRNA_cell_major_lineage_scores.tsv"), sep = "\t", na = "NA")
props <- rbindlist(sample_rows, fill = TRUE)
fwrite(props, file.path(out_dir, "PDAC_scRNA_patient_lineage_summary_long.tsv"), sep = "\t", na = "NA")

denom <- cells[, .(total_singlet_cells = .N), by = .(dataset, gsm, tissue_context, patient_id)]
wide_counts <- dcast(cells, dataset + gsm + tissue_context + patient_id ~ major_lineage, fun.aggregate = length, value.var = "cell", fill = 0)
wide <- merge(denom, wide_counts, by = c("dataset", "gsm", "tissue_context", "patient_id"), all = TRUE)
lineage_cols <- setdiff(names(wide), c("dataset", "gsm", "tissue_context", "patient_id", "total_singlet_cells"))
for (nm in lineage_cols) wide[[paste0("prop_", nm)]] <- wide[[nm]] / wide$total_singlet_cells
module_sample <- cells[, .(malignant_epithelial_program = mean(module_malignant_epithelial_program),
                           CAF_program = mean(module_CAF_program), myeloid_program = mean(module_myeloid_program),
                           CD8_activation = mean(module_CD8_activation), CD8_exhaustion = mean(module_CD8_exhaustion)),
                       by = .(dataset, gsm, tissue_context, patient_id)]
wide <- merge(wide, module_sample, by = c("dataset", "gsm", "tissue_context", "patient_id"), all = TRUE)
fwrite(wide, file.path(out_dir, "PDAC_scRNA_patient_lineage_summary.tsv"), sep = "\t", na = "NA")

primary <- wide[tissue_context == "primary_tumor"]
cor_rows <- list()
for (ds in unique(primary$dataset)) {
  x <- primary[dataset == ds]
  if (nrow(x) < 6) next
  safe_cor <- function(a, b) if (length(unique(a)) > 2 && length(unique(b)) > 2) suppressWarnings(cor(a, b, method = "spearman", use = "pairwise.complete.obs")) else NA_real_
  cor_rows <<- append(cor_rows, list(data.frame(dataset = ds, n_patients = nrow(x),
    rho_malignant_CAF = safe_cor(x$malignant_epithelial_program, x$CAF_program),
    rho_malignant_myeloid = safe_cor(x$malignant_epithelial_program, x$myeloid_program),
    rho_malignant_CD8_exhaustion = safe_cor(x$malignant_epithelial_program, x$CD8_exhaustion),
    rho_CAF_CD8_activation = safe_cor(x$CAF_program, x$CD8_activation), stringsAsFactors = FALSE)))
}
if (length(cor_rows)) fwrite(rbindlist(cor_rows), file.path(out_dir, "PDAC_scRNA_primary_ecology_correlations.tsv"), sep = "\t", na = "NA")
if (length(errors)) fwrite(rbindlist(errors, fill = TRUE), file.path(out_dir, "PDAC_scRNA_annotation_errors.tsv"), sep = "\t", na = "NA")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "annotation_sessionInfo.txt"))
