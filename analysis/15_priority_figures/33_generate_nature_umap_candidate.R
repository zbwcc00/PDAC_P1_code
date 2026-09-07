base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures/scrna_additional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(ggrepel); library(patchwork)})
suppressPackageStartupMessages({library(svglite); library(ragg)})

theme_nature <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(axis.title = element_text(size = 10), axis.text = element_text(colour = "#222222"),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, colour = "#555555", hjust = 0.5),
        legend.position = "none", panel.border = element_rect(colour = "#222222", fill = NA, linewidth = 0.35),
        plot.margin = margin(8, 8, 8, 8))
pal <- c(endothelial = "#0072B2", fibroblast_CAF = "#56B4E9", myeloid = "#009E73",
         T_NK = "#CC79A7", B_cell = "#E69F00", epithelial = "#999999",
         malignant_epithelial = "#264653", acinar = "#F0E442", mast = "#D55E00",
         unknown_ambiguous = "#BDBDBD")
dataset_files <- c(GSE154778 = "GSE154778__primary_tumor.rds", GSE155698 = "GSE155698__primary_tumor.rds", GSE212966 = "GSE212966__primary_tumor.rds")

make_umap <- function(ds, file, n_show = 6000) {
  obj <- readRDS(file.path(base_dir, "analysis/03_scrna_annotation/seurat_annotated", file))
  set.seed(20260905 + match(ds, names(dataset_files)))
  keep <- sample(colnames(obj), min(n_show, ncol(obj)))
  obj <- subset(obj, cells = keep)
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
  lab <- merge(lab, counts, by = "lineage"); lab$label <- paste0(as.character(lab$lineage), "\n(n=", format(lab$n, big.mark = ","), ")")
  lab$label <- gsub("_", " ", lab$label, fixed = TRUE)
  ggplot(emb, aes(UMAP1, UMAP2, colour = lineage)) +
    geom_point(size = 0.25, alpha = 0.42, stroke = 0) +
    ggrepel::geom_label_repel(data = lab, aes(label = label), inherit.aes = FALSE, x = lab$UMAP1, y = lab$UMAP2,
                              colour = "#222222", fill = scales::alpha("white", 0.82), label.size = 0.15,
                              size = 2.85, fontface = "bold", box.padding = 0.7, point.padding = 0.35,
                              min.segment.length = 0, seed = 20260905, max.overlaps = Inf) +
    scale_colour_manual(values = pal, drop = FALSE) + coord_equal() +
    labs(title = ds, subtitle = paste0("Primary tumour cells shown: n=", format(nrow(emb), big.mark = ",")), x = "UMAP 1", y = "UMAP 2") + theme_nature
}

p <- make_umap("GSE154778", dataset_files[["GSE154778"]]) | make_umap("GSE155698", dataset_files[["GSE155698"]]) | make_umap("GSE212966", dataset_files[["GSE212966"]])
p <- p + plot_annotation(title = "Single-cell landscapes across independent PDAC discovery cohorts", subtitle = "Directly labelled major lineages; each cohort is embedded separately and shown with a fixed random sample.") & theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5), plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "#555555"))
stem <- file.path(out_dir, "Figure_candidate_nature_style_scRNA_UMAP")
ggsave(paste0(stem, ".pdf"), p, width = 13.5, height = 4.8, device = cairo_pdf)
svglite::svglite(paste0(stem, ".svg"), width = 13.5, height = 4.8); print(p); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = 13.5, height = 4.8, units = "in", res = 600, compression = "lzw"); print(p); dev.off()
png(paste0(stem, ".png"), width = 4050, height = 1440, res = 300, type = "cairo-png"); print(p); dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "33_R_sessionInfo.txt"))
