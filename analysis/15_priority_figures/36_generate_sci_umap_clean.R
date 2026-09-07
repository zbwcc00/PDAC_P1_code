base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures/scrna_additional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(ggrastr)
  library(patchwork)
  library(svglite)
  library(ragg)
})

theme_sci <- theme_classic(base_size = 10.5, base_family = "Arial") +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(colour = "#222222", size = 8.5),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8.5, colour = "#555555", hjust = 0.5),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8),
        legend.key.width = unit(0.8, "cm"), panel.border = element_rect(colour = "#222222", fill = NA, linewidth = 0.35),
        plot.margin = margin(5, 6, 5, 6))

pal <- c(endothelial = "#007C91", fibroblast_CAF = "#C65A11", myeloid = "#5A4A9F",
         T_NK = "#B2185B", B_cell = "#2F7F4F", epithelial = "#A66A00",
         malignant_epithelial = "#8C3B2A", acinar = "#1769AA", mast = "#B51D1D",
         unknown_ambiguous = "#9E9E9E")
display_names <- c(endothelial = "Endothelial", fibroblast_CAF = "CAF", myeloid = "Myeloid", T_NK = "T/NK",
                   B_cell = "B cell", epithelial = "Epithelial", malignant_epithelial = "Malignant epithelial",
                   acinar = "Acinar", mast = "Mast", unknown_ambiguous = "Unresolved")
dataset_files <- c(GSE154778 = "GSE154778__primary_tumor.rds", GSE155698 = "GSE155698__primary_tumor.rds", GSE212966 = "GSE212966__primary_tumor.rds")

make_umap <- function(ds, file, n_show = 6000) {
  obj <- readRDS(file.path(base_dir, "analysis/03_scrna_annotation/seurat_annotated", file))
  set.seed(20260905 + match(ds, names(dataset_files)))
  keep <- sample(colnames(obj), min(n_show, ncol(obj))); obj <- subset(obj, cells = keep)
  if (length(VariableFeatures(obj)) == 0) obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 25, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, n.neighbors = 35, min.dist = 0.38, spread = 1.15, seed.use = 20260905, verbose = FALSE)
  emb <- as.data.frame(Embeddings(obj, "umap")); colnames(emb)[1:2] <- c("UMAP1", "UMAP2")
  emb$lineage <- as.character(obj$major_lineage); emb <- emb[!is.na(emb$lineage) & emb$lineage %in% names(pal), ]
  emb$lineage <- factor(emb$lineage, levels = names(pal))
  lab <- aggregate(cbind(UMAP1, UMAP2) ~ lineage, emb, median)
  lab <- lab[lab$lineage != "unknown_ambiguous", ]
  counts <- as.data.frame(table(emb$lineage)); names(counts) <- c("lineage", "n")
  lab <- merge(lab, counts, by = "lineage")
  lab$label <- display_names[as.character(lab$lineage)]
  ggplot(emb, aes(UMAP1, UMAP2, colour = lineage)) +
    ggrastr::geom_point_rast(size = 0.43, alpha = 0.92, stroke = 0, raster.dpi = 600) +
    ggrepel::geom_label_repel(data = lab, aes(label = label), inherit.aes = FALSE,
                              x = lab$UMAP1, y = lab$UMAP2, colour = "#161616", fill = "white",
                              alpha = 0.88, size = 2.65, fontface = "bold", label.size = 0.15,
                              label.padding = unit(0.12, "lines"), box.padding = 0.55, point.padding = 0.25,
                              segment.colour = NA, force = 2.2, max.time = 4, seed = 20260905,
                              max.overlaps = Inf) +
    scale_colour_manual(values = pal, labels = display_names, drop = FALSE) + coord_equal() +
    labs(title = ds, subtitle = paste0("Primary tumour cells shown: n=", format(nrow(emb), big.mark = ",")), x = "UMAP 1", y = "UMAP 2") +
    theme_sci
}

p <- make_umap("GSE154778", dataset_files[["GSE154778"]]) |
     make_umap("GSE155698", dataset_files[["GSE155698"]]) |
     make_umap("GSE212966", dataset_files[["GSE212966"]])
p <- p + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p <- p + plot_annotation(title = "Single-cell landscapes across independent PDAC discovery cohorts",
                         subtitle = "Independent embeddings; saturated lineage colours and compact labels without guide lines.") &
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "#555555"))

stem <- file.path(out_dir, "Figure_candidate_sci_style_scRNA_UMAP_clean")
ggsave(paste0(stem, ".pdf"), p, width = 13.5, height = 5.15, device = cairo_pdf)
svglite::svglite(paste0(stem, ".svg"), width = 13.5, height = 5.15); print(p); dev.off()
tiff(paste0(stem, ".tiff"), width = 13.5, height = 5.15, units = "in", res = 300, compression = "lzw", type = "cairo"); print(p); dev.off()
png(paste0(stem, ".png"), width = 4050, height = 1545, res = 300, type = "cairo-png"); print(p); dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "36_R_sessionInfo.txt"))
cat("Wrote", stem, "in PDF/SVG/TIFF/PNG formats\n")
