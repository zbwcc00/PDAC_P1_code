options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(png)
  library(grid)
  library(patchwork)
  library(jsonlite)
})

root <- "D:/PDAC_P1"
temp_out <- file.path(tempdir(), "pdac_spatial_figures")
out <- file.path(root, "analysis/09_spatial_validation/figures")
dir.create(temp_out, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, stem, width, height) {
  temporary_png <- file.path(temp_out, paste0(stem, ".png"))
  temporary_pdf <- file.path(temp_out, paste0(stem, ".pdf"))
  ggsave(temporary_png, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(temporary_pdf, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  file.copy(temporary_png, file.path(out, paste0(stem, ".png")), overwrite = TRUE)
  file.copy(temporary_pdf, file.path(out, paste0(stem, ".pdf")), overwrite = TRUE)
}

make_overlay <- function(score_file, image_file, scale_file, sample_id, feature, title_text) {
  scores <- fread(score_file)[sample == sample_id & in_tissue == 1]
  scale_info <- fromJSON(gzfile(scale_file))
  scale_value <- as.numeric(scale_info$tissue_hires_scalef)
  image <- readPNG(image_file)
  image_height <- dim(image)[1]
  image_width <- dim(image)[2]
  scores[, `:=`(x = pxl_col_in_fullres * scale_value, y = image_height - pxl_row_in_fullres * scale_value, value = get(feature))]
  scores <- scores[is.finite(value)]
  ggplot() + annotation_raster(image, xmin = 0, xmax = image_width, ymin = 0, ymax = image_height, interpolate = TRUE) + geom_point(data = scores, aes(x = x, y = y, color = value), size = 0.22, alpha = 0.72) + scale_color_viridis_c(option = "magma", name = feature) + coord_fixed(xlim = c(0, image_width), ylim = c(0, image_height), expand = FALSE) + theme_void() + theme(legend.position = "right", plot.title = element_text(size = 9, face = "bold")) + ggtitle(title_text)
}

make_gse282_overlay <- function(sample_id, feature) {
  raw <- file.path(root, "data/02_scrna_external/GSE282302/raw")
  score_file <- file.path(root, "analysis/09_spatial_validation/GSE282302/GSE282302_spot_anchor_scores.tsv")
  image_file <- file.path(root, "analysis/09_spatial_validation/GSE282302/tissue_images_uncompressed", paste0(sample_id, "_tissue_hires_image.png"))
  scale_file <- file.path(raw, paste0(sample_id, "_scalefactors_json.json.gz"))
  patient_id <- sub(".*_(D[0-9]+)_.*$", "\\1", sample_id)
  make_overlay(score_file, image_file, scale_file, sample_id, feature, paste0(patient_id, " | ", feature))
}

roi <- fread(file.path(root, "analysis/09_spatial_validation/GSE282302/GSE282302_ROI_spatial_anchor_results.tsv"))
selected <- roi[order(patient, sample), .SD[1], by = patient]$sample
selected <- selected[seq_len(min(14, length(selected)))]
features <- c("EPAS1_logCPM", "vascular_marker_score", "endothelial_program_excluding_EPAS1")
for (sample_id in selected) {
  panels <- lapply(features, function(feature) make_gse282_overlay(sample_id, feature))
  save_plot(wrap_plots(plotlist = panels, ncol = 3), paste0("GSE282302_overlay_", sample_id), 12, 4)
}

patient <- fread(file.path(root, "analysis/09_spatial_validation/GSE282302/GSE282302_patient_spatial_anchor_results.tsv"))
patient[, patient := factor(patient, levels = patient[order(median_EPAS1_vascular_rho)])]
correlation_plot <- ggplot(patient, aes(x = median_EPAS1_vascular_rho, y = patient)) + geom_segment(aes(x = 0, xend = median_EPAS1_vascular_rho, y = patient, yend = patient), color = "grey70", linewidth = 0.7) + geom_point(size = 2.4, color = "#0072B2") + geom_vline(xintercept = 0, linetype = 2, color = "grey40") + theme_classic(base_size = 10) + labs(x = "Patient-median Spearman rho: EPAS1 vs vascular markers", y = NULL) + theme(axis.text.y = element_text(size = 8))
program_plot <- ggplot(patient, aes(x = median_EPAS1_program_rho, y = patient)) + geom_segment(aes(x = 0, xend = median_EPAS1_program_rho, y = patient, yend = patient), color = "grey70", linewidth = 0.7) + geom_point(size = 2.4, color = "#D55E00", na.rm = TRUE) + geom_vline(xintercept = 0, linetype = 2, color = "grey40") + theme_classic(base_size = 10) + labs(x = "Patient-median rho: EPAS1 vs endothelial programme", y = NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save_plot(correlation_plot + program_plot + plot_layout(widths = c(1.2, 1)), "GSE282302_patient_level_spatial_correlations", 9, 5)

representative <- selected[seq_len(min(3, length(selected)))]
main_panels <- unlist(lapply(representative, function(sample_id) lapply(c("EPAS1_logCPM", "vascular_marker_score"), function(feature) make_gse282_overlay(sample_id, feature))), recursive = FALSE)
main_figure <- wrap_plots(main_panels, ncol = 3) / (correlation_plot + program_plot + plot_layout(widths = c(1.2, 1))) + plot_annotation(title = "Spatial anchoring of EPAS1 to vascular-marker regions in PDAC", subtitle = "Independent multi-patient Visium validation; dots are spots overlaid on the tissue image")
save_plot(main_figure, "Figure_main_spatial_EPAS1_vascular_anchor", 13, 9)

writeLines(c(
  "# Spatial image-overlay figure manifest",
  "",
  paste0("- GSE282302 representative patient panels: ", length(selected), "."),
  "- Every overlay uses the original tissue_hires_image, tissue_hires_scalef, and full-resolution spot coordinates.",
  "- EPAS1/vascular-marker overlays are the primary figure; the endothelial programme excluding EPAS1 is shown as a sensitivity panel.",
  "- Main figure: `Figure_main_spatial_EPAS1_vascular_anchor.pdf` and `.png`.",
  "- Interpretation: spatial EPAS1-to-vascular anchoring is replicated; the complete differential programme is not uniformly reproduced."
), file.path(out, "spatial_image_overlay_figure_manifest.md"))
