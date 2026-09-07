root <- "D:/PDAC_P1"
fig_dir <- file.path(root, "analysis/15_priority_figures")
out_dir <- file.path(fig_dir, "main_figures_relayout")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
  library(svglite)
  library(ragg)
})

umap_file <- file.path(fig_dir, "scrna_additional/Figure_candidate_sci_style_scRNA_UMAP_clean.png")
stopifnot(file.exists(umap_file))

# Remove only the global title from the raster UMAP; keep cohort titles and legend.
umap <- png::readPNG(umap_file)
umap <- umap[180:nrow(umap), , , drop = FALSE]
umap_grob <- grid::rasterGrob(umap, interpolate = TRUE)
panel_a <- patchwork::wrap_elements(full = umap_grob)

de_files <- c(
  "GSE155698 edgeR" = file.path(root, "analysis/06_pseudobulk_DE/edgeR/GSE155698/endothelial_DE.tsv"),
  "GSE155698 DESeq2" = file.path(root, "analysis/06_pseudobulk_DE/DESeq2/GSE155698/endothelial_DE.tsv"),
  "GSE212966 edgeR" = file.path(root, "analysis/06_pseudobulk_DE/edgeR/GSE212966/endothelial_DE.tsv"),
  "GSE212966 DESeq2" = file.path(root, "analysis/06_pseudobulk_DE/DESeq2/GSE212966/endothelial_DE.tsv")
)
stopifnot(all(file.exists(de_files)))
candidates <- c("EPAS1", "TACC1", "MARCKS", "HERPUD1")

read_candidate <- function(path, label) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1] <- "gene"
  d$gene <- gsub('^"|"$', "", d$gene)
  fc_col <- if ("logFC" %in% names(d)) "logFC" else "log2FoldChange"
  fdr_col <- if ("FDR" %in% names(d)) "FDR" else "padj"
  d$logFC <- as.numeric(d[[fc_col]])
  d$FDR <- as.numeric(d[[fdr_col]])
  d <- d[d$gene %in% candidates, c("gene", "logFC", "FDR"), drop = FALSE]
  d$dataset_method <- label
  d
}
de <- do.call(rbind, Map(read_candidate, unname(de_files), names(de_files)))
grid_df <- expand.grid(gene = candidates, dataset_method = names(de_files), stringsAsFactors = FALSE)
grid_df <- merge(grid_df, de[, c("gene", "dataset_method", "logFC", "FDR")],
                 by = c("gene", "dataset_method"), all.x = TRUE, sort = FALSE)
grid_df$gene <- factor(grid_df$gene, levels = rev(candidates))
grid_df$dataset_method <- factor(grid_df$dataset_method, levels = names(de_files))
grid_df$label <- ifelse(is.na(grid_df$logFC), "NA", sprintf("%+.2f", grid_df$logFC))

lim <- max(abs(grid_df$logFC), na.rm = TRUE)
lim <- max(2, ceiling(lim * 2) / 2)
panel_b <- ggplot(grid_df, aes(dataset_method, gene, fill = logFC)) +
  geom_tile(colour = "white", linewidth = 0.65, na.rm = FALSE) +
  geom_text(aes(label = label), colour = "#111111", size = 3.35, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), na.value = "#E5E5E5",
                       name = "log2 fold-change") +
  scale_x_discrete(labels = c("GSE155698\nedgeR", "GSE155698\nDESeq2",
                              "GSE212966\nedgeR", "GSE212966\nDESeq2")) +
  labs(x = NULL, y = NULL,
       title = "Cross-cohort and cross-method convergence of endothelial candidates",
       subtitle = "Values are endothelial primary tumour versus adjacent normal log2 fold-changes") +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, colour = "#555555", hjust = 0.5),
        axis.text.x = element_text(size = 9, colour = "#222222", lineheight = 0.9),
        axis.text.y = element_text(size = 10.5, face = "bold", colour = "#222222"),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "#222222", fill = NA, linewidth = 0.35),
        legend.position = "right", legend.title = element_text(size = 9),
        legend.text = element_text(size = 8), plot.margin = margin(8, 12, 8, 12))

fig <- (panel_a / panel_b) +
  plot_layout(heights = c(1.35, 0.85), guides = "keep") +
  plot_annotation(
    title = "Figure 1. Multi-cohort single-cell discovery of the EPAS1 endothelial axis",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                  plot.tag = element_text(face = "bold", size = 15),
                  plot.margin = margin(8, 8, 8, 8))
  )

stem <- file.path(out_dir, "Figure_1_scRNA_discovery_updated")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 13.5, height = 10.0, family = "Arial")
print(fig)
dev.off()
svglite::svglite(paste0(stem, ".svg"), width = 13.5, height = 10.0); print(fig); dev.off()
tiff(paste0(stem, ".tiff"), width = 13.5, height = 10.0, units = "in", res = 600, compression = "lzw", type = "cairo"); print(fig); dev.off()
png(paste0(stem, ".png"), width = 13.5, height = 10.0, units = "in", res = 300, type = "cairo-png"); print(fig); dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "35_R_sessionInfo.txt"))
cat("Wrote", stem, "in PDF/SVG/TIFF/PNG formats\n")
