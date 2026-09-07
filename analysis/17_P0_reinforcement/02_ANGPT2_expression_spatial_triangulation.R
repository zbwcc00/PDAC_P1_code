options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})

root <- "D:/PDAC_P1"
outdir <- file.path(root, "analysis/17_P0_reinforcement")
objdir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")
datasets <- c(
  GSE154778 = "GSE154778__primary_tumor.rds",
  GSE155698 = "GSE155698__primary_tumor.rds",
  GSE212966 = "GSE212966__primary_tumor.rds"
)
genes <- c("ANGPT2", "TEK", "ITGA5", "ITGB1")
receiver_types <- c("myeloid", "T_NK", "CD8", "B_cell")

expression_rows <- list()
for (dataset in names(datasets)) {
  object <- readRDS(file.path(objdir, datasets[[dataset]]))
  metadata <- object@meta.data
  assay <- if ("RNA" %in% names(object@assays)) "RNA" else DefaultAssay(object)
  log_data <- GetAssayData(object, assay = assay, layer = "data")
  available_genes <- intersect(genes, rownames(log_data))
  for (cell_type in unique(as.character(metadata$major_lineage))) {
    cells_type <- rownames(metadata)[metadata$major_lineage == cell_type]
    if (!length(cells_type)) next
    patient_groups <- split(cells_type, as.character(metadata[cells_type, "patient_id"]))
    for (patient in names(patient_groups)) {
      cells <- patient_groups[[patient]]
      if (length(cells) < 5) next
      for (gene in available_genes) {
        values <- as.numeric(log_data[gene, cells])
        expression_rows[[length(expression_rows) + 1L]] <- data.table(
          dataset = dataset,
          patient = patient,
          cell_type = cell_type,
          gene = gene,
          n_cells = length(cells),
          mean_log_expression = mean(values),
          detection_fraction = mean(values > 0)
        )
      }
    }
  }
}
expression <- rbindlist(expression_rows)
fwrite(expression, file.path(outdir, "ANGPT2_receptor_patient_celltype_expression.tsv"), sep = "\t")

# Pre-specified evidence predicates. This is expression/proximity triangulation,
# not an independent inference that establishes a ligand-receptor interaction.
summary_rows <- list(
  expression[gene == "ANGPT2" & cell_type == "endothelial",
    .(n_patient_groups = .N, median_mean_expression = median(mean_log_expression),
      median_detection_fraction = median(detection_fraction),
      positive_patient_groups = sum(detection_fraction > 0)),
    by = .(dataset, gene, cell_type)],
  expression[gene %in% c("ITGA5", "ITGB1", "TEK") & cell_type %in% receiver_types,
    .(n_patient_groups = .N, median_mean_expression = median(mean_log_expression),
      median_detection_fraction = median(detection_fraction),
      positive_patient_groups = sum(detection_fraction > 0)),
    by = .(dataset, gene, cell_type)]
)
fwrite(rbindlist(summary_rows, fill = TRUE), file.path(outdir, "ANGPT2_receptor_expression_summary.tsv"), sep = "\t")

# Spatial evidence was generated earlier using 1,000 within-ROI permutations.
spatial_summary <- fread(file.path(root, "analysis/09_spatial_validation/SecAct/SecAct_patient_level_summary.tsv"))
angpt2_spatial <- spatial_summary[secreted == "ANGPT2"]
fwrite(angpt2_spatial, file.path(outdir, "ANGPT2_SecAct_spatial_patient_replication.tsv"), sep = "\t")

cellchat <- fread(file.path(root, "analysis/09_spatial_validation/SecAct/SecAct_candidate_CellChat_summary.tsv"))
angpt2_cellchat <- cellchat[ligand == "ANGPT2"]
fwrite(angpt2_cellchat, file.path(outdir, "ANGPT2_CellChat_candidate_edges.tsv"), sep = "\t")

report <- c(
  "# ANGPT2 communication triangulation",
  "",
  "This package adds orthogonal expression and spatial-context checks to the existing CellChat candidate edges.",
  "It does not establish direct signaling, receptor occupancy, or causality.",
  "",
  "- Expression is summarized per patient x cell type; individual cells are not independent replicates.",
  "- ANGPT2 source expression is assessed in endothelial cells.",
  "- ITGA5, ITGB1, and TEK are assessed in pre-specified immune receiver populations.",
  "- Spatial support is drawn from the prior 1,000-permutation SecAct patient-level analysis.",
  "- The final language remains: candidate endothelial-immune interface."
)
writeLines(report, file.path(outdir, "ANGPT2_triangulation_readme.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "02_sessionInfo.txt"))

