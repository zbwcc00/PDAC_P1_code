base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(png); library(grid)})
theme_pub <- theme_classic(base_size = 9, base_family = "Arial") + theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8.5, colour = "#555555"), axis.title = element_text(size = 9), axis.text = element_text(colour = "#222222"), legend.position = "bottom", legend.title = element_blank(), strip.background = element_rect(fill = "#F2F4F6", colour = NA), strip.text = element_text(face = "bold"))
save_pub <- function(p, stem, width, height) {ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf); png(file.path(out_dir, paste0(stem, ".png")), width = round(width * 300), height = round(height * 300), res = 300, type = "cairo-png"); print(p); dev.off()}
read_tsv <- function(x) read.delim(x, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")

immune <- read_tsv(file.path(base_dir, "analysis/12_functional_immune/TCGA_program_immune_GSVA_correlations.tsv"))
immune$program <- factor(immune$program, levels = c("endothelial", "T_NK")); immune$pathway <- factor(immune$pathway, levels = rev(unique(immune$pathway)))
p19a <- ggplot(immune, aes(program, pathway, fill = rho)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = 2.8) + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-0.8, 0.8)) + labs(title = "Immune-program correlations", x = NULL, y = NULL, fill = "Spearman rho") + theme_pub

gsea <- read_tsv(file.path(base_dir, "analysis/12_functional_immune/TCGA_program_correlated_Hallmark_GSEA.tsv"))
gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]; gsea <- gsea[seq_len(min(12, nrow(gsea))), ]; gsea$program <- factor(gsea$program)
gsea$term <- paste(gsea$program, gsea$pathway, sep = ": ")
gsea$term <- factor(gsea$term, levels = rev(unique(gsea$term)))
p19b <- ggplot(gsea, aes(NES, term, colour = program, size = -log10(padj))) + geom_point() + geom_vline(xintercept = 0, linetype = 2, colour = "#888888") + scale_colour_manual(values = c("T_NK" = "#CC79A7", "endothelial" = "#0072B2")) + labs(title = "Hallmark GSEA", x = "Normalized enrichment score", y = NULL, size = "-log10 adjusted P") + theme_pub

wgcna <- read_tsv(file.path(base_dir, "analysis/12_functional_immune/WGCNA/TCGA_WGCNA_module_trait_correlations.tsv"))
wgcna$module <- factor(wgcna$module, levels = unique(wgcna$module)); wgcna$trait <- factor(wgcna$trait, levels = unique(wgcna$trait))
p19c <- ggplot(wgcna, aes(trait, module, fill = rho)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = 2.4) + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) + labs(title = "WGCNA module-trait correlations", x = NULL, y = NULL, fill = "Spearman rho") + theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
p19 <- (p19a | p19b) / p19c + plot_annotation(title = "Supplementary Figure S19. Functional immune and WGCNA support", subtitle = "TCGA program-level associations support a vascular-CAF-myeloid-T/NK microenvironment; WGCNA remains associative.") & theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.5))
save_pub(p19, "Figure_S19_functional_immune_WGCNA", 8.4, 7.0)

qc <- read_tsv(file.path(base_dir, "analysis/01_tcga_qc/TCGA-PAAD_primary_tumor_qc.tsv")); qc$library_million <- qc$library_size / 1e6; qc$genes_thousand <- qc$genes_detected / 1e3
p20a <- ggplot(qc, aes(x = "Library size", y = library_million)) + geom_boxplot(fill = "#0072B2", alpha = 0.7, width = 0.45) + geom_jitter(width = 0.12, size = 1.1, alpha = 0.35) + labs(title = "Library-size QC", x = NULL, y = "Reads (million)") + theme_pub
p20b <- ggplot(qc, aes(x = "Genes detected", y = genes_thousand)) + geom_boxplot(fill = "#009E73", alpha = 0.7, width = 0.45) + geom_jitter(width = 0.12, size = 1.1, alpha = 0.35) + labs(title = "Gene-detection QC", x = NULL, y = "Genes (thousand)") + theme_pub
pca <- read_tsv(file.path(base_dir, "analysis/01_tcga_qc/TCGA-PAAD_PCA_coordinates.tsv")); pca$sample_id <- sub("-01A$", "", pca$sample_id)
p20c <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = 1.8, colour = "#D55E00", alpha = 0.7) + labs(title = "TCGA expression PCA", x = "PC1", y = "PC2") + theme_pub
p20 <- (p20a | p20b) / p20c + plot_annotation(title = "Supplementary Figure S20. TCGA-PAAD bulk quality-control overview", subtitle = "The bulk cohort contains 178 primary tumors; distributions and PCA document downstream program-scoring input quality.") & theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.5))
save_pub(p20, "Figure_S20_TCGA_bulk_QC_PCA", 7.6, 6.0)

image_files <- list.files(file.path(base_dir, "analysis/09_spatial_validation/figures"), pattern = "GSE282302_overlay_.*\\.png$", full.names = TRUE)
roi_num <- suppressWarnings(as.integer(sub(".*_D([0-9]+)_.*", "\\1", basename(image_files))))
image_files <- image_files[order(roi_num, na.last = TRUE)]
image_files <- image_files[seq_len(min(6, length(image_files)))]
draw_montage <- function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c(0.08, 0.30, 0.30, 0.30), "npc"))))
  grid.text("Supplementary Figure S21. Image-resolved EPAS1/vascular-marker spatial montage", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontfamily = "sans", fontsize = 15, fontface = "bold"))
  for (i in seq_along(image_files)) {
    img <- readPNG(image_files[[i]]); row <- ceiling(i / 2) + 1; col <- ifelse(i %% 2 == 1, 1, 2)
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col, clip = "on")); grid.raster(img, interpolate = TRUE); grid.text(sub("GSE282302_overlay_", "", basename(image_files[[i]]), fixed = TRUE), x = unit(0.02, "npc"), y = unit(0.02, "npc"), just = c("left", "bottom"), gp = gpar(fontfamily = "sans", fontsize = 6, col = "white")); popViewport()
  }
  popViewport()
}
pdf(file.path(out_dir, "Figure_S21_spatial_image_montage.pdf"), width = 8.5, height = 7.0); draw_montage(); dev.off()
png(file.path(out_dir, "Figure_S21_spatial_image_montage.png"), width = 2550, height = 2100, res = 300, type = "cairo-png"); draw_montage(); dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "26_R_sessionInfo.txt"))
