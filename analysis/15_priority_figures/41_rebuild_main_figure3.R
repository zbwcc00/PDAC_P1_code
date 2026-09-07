options(stringsAsFactors = FALSE)
root <- "D:/PDAC_P1"
out_dir <- file.path(root, "analysis/15_priority_figures/main_figures_unified")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(png); library(grid); library(jsonlite); library(svglite)
})
navy <- "#1F3B5B"; teal <- "#138A8A"; coral <- "#C85C5C"; gold <- "#C58A24"
theme_pub <- theme_classic(base_family = "Arial", base_size = 8) +
  theme(axis.line = element_line(linewidth = .35, colour = "black"), axis.ticks = element_line(linewidth = .35, colour = "black"), axis.text = element_text(colour = "black"), plot.title = element_text(face = "bold", size = 9, colour = "#111111"), plot.subtitle = element_text(size = 7, colour = "#444444"), legend.title = element_text(size = 7), legend.text = element_text(size = 6), panel.grid = element_blank())
save_pub <- function(p, stem, width = 13.5, height = 9) {
  png_path <- paste0(stem, ".png")
  png(png_path, width = width, height = height, units = "in", res = 600, type = "cairo-png"); print(p); dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "Arial"); print(p); dev.off()
  tiff(paste0(stem, ".tiff"), width = width, height = height, units = "in", res = 300, compression = "lzw", type = "cairo"); print(p); dev.off()
  svglite(paste0(stem, ".svg"), width = width, height = height); print(p); dev.off()
}

sec <- fread(file.path(root, "analysis/09_spatial_validation/SecAct/SecAct_five_targets_patient_summary_1000.tsv"))
targets <- c("SPARCL1", "VWF", "MMRN2", "ANGPT2", "IL33", "IGFBP7")
sec <- sec[secreted %in% targets]; sec$secreted <- factor(sec$secreted, levels = rev(targets)); sec$dataset <- factor(sec$dataset, levels = c("GSE282302", "GSE297144"))
pA <- ggplot(sec, aes(dataset, secreted, fill = positive_fraction)) + geom_tile(colour = "white", linewidth = .5) + geom_text(aes(label = sprintf("%d/%d", positive_patients, n_patients)), size = 2.7, fontface = "bold", colour = "#111111") + scale_fill_gradientn(colours = c("#EEF2F6", "#8AB5C7", navy), limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", ".5", "1")) + labs(title = "SecAct direction replication", subtitle = "Patient/sample-level positive direction; numerator/denominator shown", x = NULL, y = NULL, fill = "Positive fraction") + theme_pub + theme(axis.text.x = element_text(face = "bold"))

edges <- fread(file.path(root, "analysis/09_spatial_validation/SecAct/SecAct_candidate_CellChat_edges.tsv")); edges <- edges[ligand == "ANGPT2"]
bub <- edges[, .(n_edges = .N, mean_prob = mean(prob), min_p = min(pval)), by = .(source, target, ligand, receptor, pathway_name)]; bub$target <- factor(bub$target, levels = c("B_cell", "CD8", "myeloid"))
pB <- ggplot(bub, aes(ligand, target, size = n_edges, colour = mean_prob)) + geom_point(alpha = .9) + scale_size_continuous(range = c(4, 11), breaks = c(1, 3), name = "Inferred edges") + scale_colour_viridis_c(option = "C", trans = "log10", name = "Mean probability") + labs(title = "ANGPT2 CellChat contribution", subtitle = "Endothelial source -> immune receivers via ITGA5-ITGB1", x = NULL, y = NULL) + theme_pub + theme(axis.text.y = element_text(face = "bold"))

sp <- fread(file.path(root, "analysis/09_spatial_validation/GSE282302/spatial_ecology/GSE282302_spatial_ecology_spots.tsv"))
sample_id <- "GSM8641019_C1_D6_ROI1_s1"
one <- sp[sample == sample_id & is.finite(myeloid_neighbor_score) & is.finite(T_NK_neighbor_score)]
img_file <- file.path(root, "analysis/09_spatial_validation/GSE282302/tissue_images_uncompressed", paste0(sample_id, "_tissue_hires_image.png"))
scale_file <- file.path(root, "data/02_scrna_external/GSE282302/raw", paste0(sample_id, "_scalefactors_json.json.gz"))
img <- readPNG(img_file); sf <- fromJSON(gzfile(scale_file)); scale_value <- as.numeric(sf$tissue_hires_scalef); ih <- dim(img)[1]; iw <- dim(img)[2]
one[, `:=`(x = pxl_col_in_fullres * scale_value, y = ih - pxl_row_in_fullres * scale_value)]
spatial_panel <- function(dat, feature, label, pal) {
  ggplot() + annotation_raster(img, xmin = 0, xmax = iw, ymin = 0, ymax = ih, interpolate = TRUE) + geom_point(data = dat, aes(x, y, colour = get(feature)), size = .24, alpha = .78) + scale_colour_gradientn(colours = pal, name = label) + coord_fixed(xlim = c(0, iw), ylim = c(0, ih), expand = FALSE) + theme_void() + ggtitle(label) + theme(plot.title = element_text(size = 8, face = "bold", hjust = .5), legend.title = element_text(size = 6), legend.text = element_text(size = 5))
}
pC1 <- spatial_panel(one, "myeloid_neighbor_score", "Myeloid neighborhood score", c("#F2F2F2", "#D8A35D", navy))
pC2 <- spatial_panel(one, "T_NK_neighbor_score", "T/NK neighborhood score", c("#F2F2F2", "#D48787", coral))
pC <- wrap_elements(full = patchwork::patchworkGrob(pC1 + pC2 + plot_layout(guides = "collect") + plot_annotation(title = paste("Spatial neighborhood context |", sample_id), subtitle = "Representative GSE282302 ROI; spot-level visualization", theme = theme(plot.title = element_text(size = 8, hjust = .5), plot.subtitle = element_text(size = 6, hjust = .5)))))

eco <- fread(file.path(root, "analysis/09_spatial_validation/GSE282302/spatial_ecology/GSE282302_spatial_ecology_patient.tsv"))
dlong <- melt(eco, id.vars = "patient", measure.vars = c("median_myeloid_delta", "median_T_NK_delta", "median_CAF_delta", "median_malignant_delta"), variable.name = "program", value.name = "delta")
dlong$program <- factor(sub("median_", "", sub("_delta", "", dlong$program)), levels = c("myeloid", "T_NK", "CAF", "malignant"))
pD <- ggplot(dlong, aes(program, patient, fill = delta)) + geom_tile(colour = "white", linewidth = .3) + geom_text(aes(label = sprintf("%.2f", delta)), size = 2.2) + scale_fill_gradient2(low = "#3B4F7C", mid = "white", high = coral, midpoint = 0) + labs(title = "Patient-level neighborhood enrichment", subtitle = "EPAS1-high minus EPAS1-low neighborhood score", x = NULL, y = "Patient", fill = "delta score") + theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"))

fig3 <- (pA | pB) / (pC | pD) + plot_annotation(title = "Figure 3. Secreted signaling and immune spatial ecology", tag_levels = "A", theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = .5)))
stem <- file.path(out_dir, "Figure_3_secact_immune_ecology_unified")
save_pub(fig3, stem, width = 13.5, height = 10.2)
writeLines(c(
  "# Figure 3 panel audit", "",
  "- A: SecAct patient/sample-level positive-fraction heatmap across GSE282302 and GSE297144.",
  "- B: ANGPT2 endothelial-source CellChat bubble plot; duplicate myeloid edges are aggregated and represented by edge count.",
  paste0("- C: Representative GSE282302 ROI ", sample_id, " showing myeloid and T/NK neighborhood scores over the tissue image."),
  "- D: Patient-level neighborhood enrichment heatmap; each row is a patient and spot-level values are not treated as independent replicates.",
  "- No activity-versus-downstream-TF panel was added because SecAct output here does not contain a validated TF activity layer.",
  "- No chord/circle/signaling-role panel was added because the retained ANGPT2 network has one ligand, one receptor complex, and three receiver classes; a bubble plot is more faithful and less visually inflationary."
), file.path(root, "analysis/15_priority_figures/Figure3_panel_audit.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "41_R_sessionInfo.txt"))
message("Figure 3 rebuilt: ", stem)
