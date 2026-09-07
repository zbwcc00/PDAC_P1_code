source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- pdac_paths()$project_root
out_dir <- file.path(root, "analysis/15_priority_figures/main_figures_unified")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(patchwork)
  library(ggplot2)
  library(svglite)
})

navy <- "#1F3B5B"
black <- "#111111"
theme_frame <- theme_void(base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 16, colour = black, hjust = 0.5),
        plot.tag = element_text(face = "bold", size = 13, colour = navy),
        plot.margin = margin(5, 8, 5, 8))

trim_raster <- function(path, crop_top = 0, threshold = 0.985, pad = 18) {
  x <- png::readPNG(path)
  if (crop_top > 0) x <- x[(crop_top + 1):dim(x)[1], , , drop = FALSE]
  rgb <- x[, , 1:3, drop = FALSE]
  ink <- apply(rgb, c(1, 2), function(v) any(v < threshold))
  rows <- which(rowSums(ink) > 0); cols <- which(colSums(ink) > 0)
  if (!length(rows) || !length(cols)) return(grid::rasterGrob(x, interpolate = TRUE))
  r1 <- max(1, min(rows) - pad); r2 <- min(dim(x)[1], max(rows) + pad)
  c1 <- max(1, min(cols) - pad); c2 <- min(dim(x)[2], max(cols) + pad)
  grid::rasterGrob(x[r1:r2, c1:c2, , drop = FALSE], interpolate = TRUE)
}

panel <- function(path, crop_top = 0) wrap_elements(full = trim_raster(path, crop_top))

save_patchwork <- function(plot, stem, width = 13.5, height = 8.5) {
  pdf <- paste0(stem, ".pdf")
  grDevices::cairo_pdf(pdf, width = width, height = height, family = "Arial")
  print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = width, height = height)
  print(plot); dev.off()
  grDevices::tiff(paste0(stem, ".tiff"), width = width, height = height, units = "in", res = 600,
                  compression = "lzw", type = "cairo")
  print(plot); dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 300,
                 type = "cairo-png")
  print(plot); dev.off()
}

make_vertical <- function(title, panels, heights, width = 13.5, height = 8.5) {
  combined <- Reduce(`/`, panels)
  combined + plot_layout(heights = heights) +
    plot_annotation(title = title, tag_levels = "A", theme = theme_frame)
}

fig_dir <- file.path(root, "analysis/15_priority_figures")
spatial_dir <- file.path(root, "analysis/09_spatial_validation/figures")
closed_dir <- file.path(root, "analysis/13_virtual_perturbation/closed_loop")

f2 <- make_vertical(
  "Figure 2. External bulk and spatial validation of the EPAS1 vascular axis",
  list(panel(file.path(fig_dir, "Figure_S8_bulk_program_validation.png"), 92),
       panel(file.path(spatial_dir, "Figure_main_spatial_EPAS1_vascular_anchor.png"), 170)),
  c(0.9, 1.55), height = 10.4)
save_patchwork(f2, file.path(out_dir, "Figure_2_bulk_spatial_validation_unified"), height = 10.4)

f3 <- make_vertical(
  "Figure 3. Secreted signaling and immune spatial ecology",
  list(panel(file.path(fig_dir, "Figure_S10_secact_angpt2_replication.png"), 88),
       panel(file.path(fig_dir, "Figure_S9_spatial_ecology_patient_level.png"), 88)),
  c(1.0, 0.9), height = 7.7)
save_patchwork(f3, file.path(out_dir, "Figure_3_secact_immune_ecology_unified"), height = 7.7)

f4 <- make_vertical(
  "Figure 4. Virtual perturbation prioritizes the endothelial EPAS1 program",
  list(panel(file.path(closed_dir, "epas1_endothelial_program_response.png"), 0)),
  1, height = 4.8)
save_patchwork(f4, file.path(out_dir, "Figure_4_virtual_perturbation_unified"), height = 4.8)

f5 <- make_vertical(
  "Figure 5. Drug reversal and structure-guided prioritization",
  list(panel(file.path(fig_dir, "Figure_S12_DrugReflector_prioritization.png"), 90),
       panel(file.path(fig_dir, "Figure_S16_EPAS1_main_docking_QC.png"), 90),
       panel(file.path(fig_dir, "Figure_S23_EPAS1_complex_binding_modes.png"), 90)),
  c(0.75, 1.05, 0.7), height = 10.2)
save_patchwork(f5, file.path(out_dir, "Figure_5_drug_structure_closure_unified"), height = 10.2)

f6_path <- file.path(closed_dir, "manuscript_figure/Figure_6_EPAS1_evidence_chain.png")
f6 <- wrap_elements(full = trim_raster(f6_path, crop_top = 650, pad = 20)) +
  plot_annotation(title = "Figure 6. EPAS1-centered closed-loop interpretation", theme = theme_frame)
save_patchwork(f6, file.path(out_dir, "Figure_6_EPAS1_closed_loop_unified"), height = 8.7)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "38_R_sessionInfo.txt"))
cat("Unified Figure 2-6 assets written to", out_dir, "\n")
