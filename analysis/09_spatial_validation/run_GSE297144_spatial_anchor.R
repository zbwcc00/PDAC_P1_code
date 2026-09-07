options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(edgeR)})
root <- "D:/PDAC_P1"
indir <- file.path(root, "data/02_scrna_external/GSE297144/raw/extracted")
outdir <- file.path(root, "analysis/09_spatial_validation/GSE297144")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

candidates <- fread(file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
endo_up <- setdiff(candidates[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id], "EPAS1")
endo_down <- setdiff(candidates[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id], "EPAS1")
vascular_genes <- c("PECAM1", "VWF", "EMCN")
score_signed <- function(log_cpm, up, down) {
  up_present <- intersect(up, rownames(log_cpm)); down_present <- intersect(down, rownames(log_cpm))
  up_score <- if (length(up_present)) colMeans(log_cpm[up_present, , drop = FALSE]) else rep(NA_real_, ncol(log_cpm))
  down_score <- if (length(down_present)) colMeans(log_cpm[down_present, , drop = FALSE]) else rep(NA_real_, ncol(log_cpm))
  if (!length(up_present) && !length(down_present)) return(rep(NA_real_, ncol(log_cpm)))
  if (!length(up_present)) return(-down_score); if (!length(down_present)) return(up_score); (up_score - down_score) / 2
}

matrix_files <- list.files(indir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
if (!length(matrix_files)) stop("No matrix files found")
rows <- vector("list", length(matrix_files)); spots <- vector("list", length(matrix_files))
for (index in seq_along(matrix_files)) {
  matrix_file <- matrix_files[index]; sample_id <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file)); message("Processing ", index, "/", length(matrix_files), ": ", sample_id)
  feature_file <- file.path(indir, paste0(sample_id, "_features.tsv.gz")); barcode_file <- file.path(indir, paste0(sample_id, "_barcodes.tsv.gz")); position_file <- file.path(indir, paste0(sub("_FRU_.*$", "", sample_id), sub("^.*FRU", "FRU", sample_id), "_tissue_positions.csv.gz"))
  position_candidates <- list.files(indir, pattern = paste0("^", sub("_.*$", "", sample_id), ".*tissue_positions\\.csv\\.gz$"), full.names = TRUE)
  if (!length(position_candidates)) position_candidates <- list.files(indir, pattern = "_tissue_positions\\.csv\\.gz$", full.names = TRUE)[index]
  position_file <- position_candidates[1]
  features <- fread(feature_file, header = FALSE); genes <- make.unique(as.character(features[[2]])); barcodes <- fread(barcode_file, header = FALSE)[[1]]
  counts <- readMM(gzfile(matrix_file)); rownames(counts) <- genes; colnames(counts) <- barcodes
  library_size <- Matrix::colSums(counts); requested <- unique(c("EPAS1", vascular_genes, endo_up, endo_down)); present <- intersect(requested, rownames(counts)); log_cpm <- log2(sweep(counts[present, , drop = FALSE], 2, library_size / 1e6, "/") + 1)
  epas1 <- if ("EPAS1" %in% rownames(log_cpm)) as.numeric(log_cpm["EPAS1", ]) else rep(NA_real_, ncol(counts)); vascular_present <- intersect(vascular_genes, rownames(log_cpm)); vascular_score <- if (length(vascular_present)) colMeans(t(scale(t(log_cpm[vascular_present, , drop = FALSE]))), na.rm = TRUE) else rep(NA_real_, ncol(counts)); endothelial_program <- score_signed(t(scale(t(log_cpm))), endo_up, endo_down)
  positions <- fread(position_file); positions <- positions[match(barcodes, barcode)]; positions[, `:=`(sample = sample_id, patient = sub(".*_(FRU_[0-9]+)_.*$", "\\1", sample_id), n_genes = Matrix::colSums(counts > 0), total_counts = library_size, EPAS1_logCPM = epas1, vascular_marker_score = vascular_score, endothelial_program_excluding_EPAS1 = endothelial_program)]
  positions[, EPAS1_quartile := fifelse(epas1 >= quantile(epas1, .75, na.rm = TRUE), "high", fifelse(epas1 <= quantile(epas1, .25, na.rm = TRUE), "low", "middle"))]
  valid <- positions[is.finite(EPAS1_logCPM) & is.finite(vascular_marker_score)]; program_valid <- positions[is.finite(EPAS1_logCPM) & is.finite(endothelial_program_excluding_EPAS1) & is.finite(vascular_marker_score)]
  gsm_prefix <- sub("_.*$", "", sample_id)
  rows[[index]] <- data.table(sample = sample_id, patient = unique(positions$patient), n_spots = ncol(counts), n_position_matched = sum(!is.na(positions$pxl_row_in_fullres)), n_genes_median = median(positions$n_genes), total_counts_median = median(positions$total_counts), n_endo_up = length(intersect(endo_up, rownames(log_cpm))), n_endo_down = length(intersect(endo_down, rownames(log_cpm))), n_vascular_markers = length(vascular_present), EPAS1_vascular_rho = if (nrow(valid) >= 20) cor(valid$EPAS1_logCPM, valid$vascular_marker_score, method = "spearman") else NA_real_, EPAS1_program_rho = if (nrow(program_valid) >= 20) cor(program_valid$EPAS1_logCPM, program_valid$endothelial_program_excluding_EPAS1, method = "spearman") else NA_real_, program_vascular_rho = if (nrow(program_valid) >= 20) cor(program_valid$endothelial_program_excluding_EPAS1, program_valid$vascular_marker_score, method = "spearman") else NA_real_, image_present = length(list.files(indir, pattern = paste0("^", gsm_prefix, ".*tissue_hires_image\\.png\\.gz$"))) > 0, scalefactor_present = file.exists(file.path(indir, paste0(sample_id, "_scalefactors_json.json.gz"))))
  spots[[index]] <- positions[, .(sample, patient, barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres, n_genes, total_counts, EPAS1_logCPM, vascular_marker_score, endothelial_program_excluding_EPAS1, EPAS1_quartile)]
}
roi <- rbindlist(rows, fill = TRUE); spot <- rbindlist(spots, fill = TRUE); patient <- roi[, .(n_ROI = .N, n_spots = sum(n_spots), median_EPAS1_vascular_rho = median(EPAS1_vascular_rho, na.rm = TRUE), median_EPAS1_program_rho = median(EPAS1_program_rho, na.rm = TRUE), proportion_positive_EPAS1_vascular = mean(EPAS1_vascular_rho > 0, na.rm = TRUE), proportion_positive_EPAS1_program = mean(EPAS1_program_rho > 0, na.rm = TRUE)), by = patient]
fwrite(roi, file.path(outdir, "GSE297144_ROI_spatial_anchor_results.tsv"), sep = "\t"); fwrite(patient, file.path(outdir, "GSE297144_patient_spatial_anchor_results.tsv"), sep = "\t"); fwrite(spot, file.path(outdir, "GSE297144_spot_anchor_scores.tsv"), sep = "\t")
writeLines(c("# GSE297144 independent spatial anchor validation", "", paste0("- Processed ", nrow(roi), " ROIs from ", uniqueN(roi$patient), " patients."), "- Matrix, spot coordinates, scalefactors, and high-resolution image availability were checked for each sample.", "- EPAS1 vascular-marker correlation is the primary spatial anchor; the differential endothelial programme excluding EPAS1 is a secondary sensitivity analysis.", "- Patient, rather than spot, is the replication unit.", "- This supports spatial molecular anchoring and does not establish causality or treatment efficacy."), file.path(outdir, "GSE297144_spatial_anchor_report.md"))
