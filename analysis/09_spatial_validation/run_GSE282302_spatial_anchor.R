options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(Seurat)
  library(ggplot2)
})

source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- pdac_paths()$project_root
indir <- file.path(root, "data/02_scrna_external/GSE282302/raw")
outdir <- file.path(root, "analysis/09_spatial_validation/GSE282302")
plotdir <- file.path(outdir, "coordinate_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

candidate_file <- file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv")
candidates <- fread(candidate_file)
endo_up <- candidates[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id]
endo_down <- candidates[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id]
endo_up <- setdiff(endo_up, "EPAS1")
endo_down <- setdiff(endo_down, "EPAS1")
vascular_genes <- c("PECAM1", "VWF", "EMCN")

score_signed <- function(log_cpm, up, down) {
  up_present <- intersect(up, rownames(log_cpm))
  down_present <- intersect(down, rownames(log_cpm))
  up_score <- if (length(up_present)) colMeans(log_cpm[up_present, , drop = FALSE]) else rep(NA_real_, ncol(log_cpm))
  down_score <- if (length(down_present)) colMeans(log_cpm[down_present, , drop = FALSE]) else rep(NA_real_, ncol(log_cpm))
  if (!length(up_present) && !length(down_present)) return(rep(NA_real_, ncol(log_cpm)))
  if (!length(up_present)) return(-down_score)
  if (!length(down_present)) return(up_score)
  (up_score - down_score) / 2
}

matrix_files <- list.files(indir, pattern = "_filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
if (!length(matrix_files)) stop("No filtered H5 matrices found in ", indir)
max_roi <- suppressWarnings(as.integer(Sys.getenv("GSE282302_MAX_ROI", "0")))
if (is.finite(max_roi) && max_roi > 0) matrix_files <- matrix_files[seq_len(min(max_roi, length(matrix_files)))]
results <- vector("list", length(matrix_files))
spot_results <- vector("list", length(matrix_files))

for (index in seq_along(matrix_files)) {
  matrix_file <- matrix_files[index]
  sample_id <- sub("_filtered_feature_bc_matrix\\.h5$", "", basename(matrix_file))
  message("Processing ", index, "/", length(matrix_files), ": ", sample_id)
  position_file <- file.path(indir, paste0(sample_id, "_tissue_positions_list.csv.gz"))
  scale_file <- file.path(indir, paste0(sample_id, "_scalefactors_json.json.gz"))
  image_file <- file.path(indir, paste0(sample_id, "_tissue_hires_image.png.gz"))
  counts <- Read10X_h5(matrix_file)
  if (is.list(counts)) counts <- counts[[1]]
  rownames(counts) <- make.unique(rownames(counts))
  library_size <- Matrix::colSums(counts)
  requested_genes <- unique(c("EPAS1", vascular_genes, endo_up, endo_down))
  requested_present <- intersect(requested_genes, rownames(counts))
  selected_counts <- counts[requested_present, , drop = FALSE]
  log_cpm <- log2(sweep(selected_counts, 2, library_size / 1e6, "/") + 1)
  epas1 <- if ("EPAS1" %in% rownames(log_cpm)) as.numeric(log_cpm["EPAS1", ]) else rep(NA_real_, ncol(log_cpm))
  vascular_present <- intersect(vascular_genes, rownames(log_cpm))
  vascular_score <- if (length(vascular_present)) colMeans(t(scale(t(log_cpm[vascular_present, , drop = FALSE]))), na.rm = TRUE) else rep(NA_real_, ncol(log_cpm))
  endothelial_program <- score_signed(t(scale(t(log_cpm))), endo_up, endo_down)
  positions <- fread(position_file, header = FALSE)
  setnames(positions, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"))
  positions <- positions[match(colnames(counts), barcode)]
  positions[, `:=`(EPAS1_logCPM = epas1, vascular_marker_score = vascular_score, endothelial_program_excluding_EPAS1 = endothelial_program)]
  positions[, `:=`(sample = sample_id, patient = sub(".*_(D[0-9]+)(?:_.*)?$", "\\1", sample_id), n_genes = Matrix::colSums(counts > 0), total_counts = Matrix::colSums(counts))]
  positions[, EPAS1_quartile := fifelse(epas1 >= quantile(epas1, 0.75, na.rm = TRUE), "high", fifelse(epas1 <= quantile(epas1, 0.25, na.rm = TRUE), "low", "middle"))]
  valid <- positions[is.finite(EPAS1_logCPM) & is.finite(vascular_marker_score)]
  program_valid <- positions[is.finite(EPAS1_logCPM) & is.finite(endothelial_program_excluding_EPAS1)]
  epas1_vascular_rho <- if (nrow(valid) >= 20) suppressWarnings(cor(valid$EPAS1_logCPM, valid$vascular_marker_score, method = "spearman")) else NA_real_
  epas1_program_rho <- if (nrow(program_valid) >= 20) suppressWarnings(cor(program_valid$EPAS1_logCPM, program_valid$endothelial_program_excluding_EPAS1, method = "spearman")) else NA_real_
  program_vascular_rho <- if (nrow(program_valid) >= 20) suppressWarnings(cor(program_valid$endothelial_program_excluding_EPAS1, program_valid$vascular_marker_score, method = "spearman")) else NA_real_
  results[[index]] <- data.table(sample = sample_id, patient = unique(positions$patient), n_spots = ncol(counts), n_in_tissue = sum(positions$in_tissue == 1, na.rm = TRUE), n_position_matched = sum(!is.na(positions$pxl_row_in_fullres)), n_genes_median = median(positions$n_genes), total_counts_median = median(positions$total_counts), n_endo_up = length(intersect(endo_up, rownames(log_cpm))), n_endo_down = length(intersect(endo_down, rownames(log_cpm))), n_vascular_markers = length(vascular_present), EPAS1_vascular_rho = epas1_vascular_rho, EPAS1_program_rho = epas1_program_rho, program_vascular_rho = program_vascular_rho, image_present = file.exists(image_file), scalefactor_present = file.exists(scale_file))
  spot_results[[index]] <- positions[, .(sample, patient, barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres, n_genes, total_counts, EPAS1_logCPM, vascular_marker_score, endothelial_program_excluding_EPAS1, EPAS1_quartile)]

  if (index <= 12) {
    plot_data <- positions[in_tissue == 1 & is.finite(EPAS1_logCPM)]
    for (feature in c("EPAS1_logCPM", "vascular_marker_score", "endothelial_program_excluding_EPAS1")) {
      plot_data[, plotted_value := get(feature)]
      plot <- ggplot(plot_data, aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres, color = plotted_value)) + geom_point(size = 0.45, alpha = 0.8) + scale_color_viridis_c(na.value = "grey80") + coord_fixed() + theme_void() + labs(title = paste(sample_id, feature), color = feature)
      temporary_plot <- file.path(tempdir(), paste0(sample_id, "_", feature, ".png"))
      ggsave(temporary_plot, plot, width = 6, height = 5, dpi = 220)
      file.copy(temporary_plot, file.path(plotdir, paste0(sample_id, "_", feature, ".png")), overwrite = TRUE)
    }
  }
}

sample_table <- rbindlist(results, fill = TRUE)
spot_table <- rbindlist(spot_results, fill = TRUE)
patient_table <- sample_table[, .(n_ROI = .N, n_spots = sum(n_spots), median_EPAS1_vascular_rho = median(EPAS1_vascular_rho, na.rm = TRUE), median_EPAS1_program_rho = median(EPAS1_program_rho, na.rm = TRUE), proportion_positive_EPAS1_vascular = mean(EPAS1_vascular_rho > 0, na.rm = TRUE), proportion_positive_EPAS1_program = mean(EPAS1_program_rho > 0, na.rm = TRUE)), by = patient]
patient_table[, `:=`(EPAS1_vascular_positive = median_EPAS1_vascular_rho > 0, EPAS1_program_positive = median_EPAS1_program_rho > 0)]

fwrite(sample_table, file.path(outdir, "GSE282302_ROI_spatial_anchor_results.tsv"), sep = "\t")
fwrite(patient_table, file.path(outdir, "GSE282302_patient_spatial_anchor_results.tsv"), sep = "\t")
fwrite(spot_table, file.path(outdir, "GSE282302_spot_anchor_scores.tsv"), sep = "\t")
writeLines(c(
  "# GSE282302 spatial anchor validation",
  "",
  paste0("- Processed ", nrow(sample_table), " spatial ROIs from ", uniqueN(sample_table$patient), " patient identifiers.") ,
  "- Each ROI was checked for filtered H5 counts, tissue-position coordinates, scalefactors, and a high-resolution tissue image.",
  "- EPAS1, PECAM1/VWF/EMCN vascular-marker score, and the candidate endothelial program excluding EPAS1 were calculated within each ROI.",
  "- Primary replication unit is the patient, not individual spots. ROI-level correlations are descriptive and are not treated as independent patient replicates.",
  "- Coordinate plots are generated for the first 12 ROIs; the original image files remain available for image-aligned visualization.",
  "- This validates spatial molecular anchoring and does not establish cell-type identity, causality, or drug efficacy."
), file.path(outdir, "GSE282302_spatial_anchor_report.md"))
