options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

root <- "D:/PDAC_P1"
out <- file.path(root, "analysis/13_virtual_perturbation/closed_loop/supplement_gene_level")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE154778", "GSE155698", "GSE212966")
result_dir <- file.path(root, "analysis/13_virtual_perturbation/results")
object_dir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")

read_one <- function(dataset) {
  result_file <- file.path(result_dir, paste0(dataset, "_primary__endothelial__EPAS1__KO.rds"))
  object_file <- file.path(object_dir, paste0(dataset, "__primary_tumor.rds"))
  res <- readRDS(result_file)
  dr <- as.data.frame(res$diffRegulation)
  dr$dataset <- dataset
  dr$padj_safe <- pmax(as.numeric(dr$p.adj), .Machine$double.xmin)
  dr$neglog10_fdr <- -log10(dr$padj_safe)
  dr$direction <- ifelse(dr$Z >= 0, "Positive", "Negative")
  dr$significant <- dr$p.adj < 0.05
  dr$rank_score <- rank(-dr$neglog10_fdr, ties.method = "min")

  obj <- readRDS(object_file)
  assay <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
  counts <- GetAssayData(obj, assay = assay, layer = "counts")
  cell_type <- as.character(obj@meta.data$major_lineage)
  cells <- rownames(obj@meta.data)[cell_type == "endothelial"]
  cells <- intersect(cells, colnames(counts))
  counts <- counts[, cells, drop = FALSE]
  lib <- Matrix::colSums(counts)
  counts <- counts[, lib >= 100, drop = FALSE]
  lib <- Matrix::colSums(counts)
  wanted <- intersect(dr$gene, rownames(counts))
  norm <- t(t(counts[wanted, , drop = FALSE]) / pmax(lib, 1)) * 1e4
  mean_expr <- Matrix::rowMeans(log1p(norm))
  expr_df <- data.frame(gene = names(mean_expr), mean_log1p_cpm = as.numeric(mean_expr))
  dr <- merge(dr, expr_df, by = "gene", all.x = TRUE, sort = FALSE)
  dr$n_endothelial_cells <- ncol(counts)
  rm(obj, counts, norm, res)
  gc(verbose = FALSE)
  dr
}

all_dr <- do.call(rbind, lapply(datasets, read_one))
all_dr$dataset <- factor(all_dr$dataset, levels = datasets)
all_dr$neglog10_fdr_plot <- pmin(all_dr$neglog10_fdr, 50)

# Cross-dataset rank table: a compact, auditable list rather than a claim of shared DE.
top_each <- do.call(rbind, lapply(datasets, function(ds) {
  z <- all_dr[all_dr$dataset == ds, , drop = FALSE]
  z <- z[order(z$rank_score, -abs(z$Z)), , drop = FALSE]
  z <- z[seq_len(min(15, nrow(z))), , drop = FALSE]
  z$within_dataset_rank <- seq_len(nrow(z))
  z
}))
write.table(top_each[, c("dataset", "within_dataset_rank", "gene", "Z", "FC", "p.value", "p.adj", "mean_log1p_cpm", "n_endothelial_cells")],
            file.path(out, "EPAS1_KO_top_genes_by_dataset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rank_summary <- do.call(rbind, lapply(split(all_dr, all_dr$gene), function(z) {
  data.frame(gene = z$gene[1], n_datasets = nrow(z), n_sig = sum(z$significant),
             median_abs_Z = median(abs(z$Z), na.rm = TRUE), median_rank = median(z$rank_score, na.rm = TRUE))
}))
rank_summary <- rank_summary[order(-rank_summary$n_sig, rank_summary$median_rank, -rank_summary$median_abs_Z), ]
write.table(rank_summary, file.path(out, "EPAS1_KO_cross_dataset_rank_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pal <- c(Positive = "#C23B4A", Negative = "#1F78A8")
sig_cols <- c(`Not significant` = "#B8C2CC", Significant = "#C23B4A")
theme_pub <- theme_classic(base_family = "Arial", base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "#EAF0F5", colour = NA),
        strip.text = element_text(face = "bold", colour = "#16324F"),
        axis.title = element_text(face = "bold"),
        legend.position = "top", panel.grid = element_blank())

# A. Volcano-like display uses signed network Z, not the unstable raw FC scale.
volc <- ggplot(all_dr, aes(x = Z, y = neglog10_fdr_plot, colour = direction, size = significant)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "#7B8794", linewidth = .35) +
  geom_vline(xintercept = 0, colour = "#7B8794", linewidth = .3) +
  geom_point(alpha = .8) +
  geom_text_repel(data = subset(all_dr, significant), aes(label = gene), family = "Arial",
                  size = 2.7, max.overlaps = 25, box.padding = .35, point.padding = .15,
                  show.legend = FALSE) +
  facet_wrap(~dataset, nrow = 1) + scale_colour_manual(values = pal) +
  scale_size_manual(values = c(`FALSE` = 1.5, `TRUE` = 2.8), guide = "none") +
  labs(title = "EPAS1 knockout gene-level network response", x = "Signed network Z",
       y = expression(-log[10](adjusted~P)~(capped~at~50)), colour = NULL) + theme_pub

# B. Top-gene rank heatmap, keeping only the 15 highest-ranked genes in any cohort.
heat_genes <- unique(c(rank_summary$gene[rank_summary$n_sig > 0],
                       unlist(lapply(datasets, function(ds) {
                         z <- top_each[top_each$dataset == ds & top_each$gene != "EPAS1", ]
                         head(z$gene, 3)
                       }))))
heat_df <- all_dr[all_dr$gene %in% heat_genes, c("dataset", "gene", "Z")]
gene_order <- unique(top_each$gene[order(top_each$within_dataset_rank)])
gene_order <- gene_order[gene_order %in% heat_genes]
heat_df$gene <- factor(heat_df$gene, levels = rev(gene_order))
rank_plot <- ggplot(heat_df, aes(x = dataset, y = gene, fill = Z)) +
  geom_tile(colour = "white", linewidth = .3) +
  geom_text(aes(label = sprintf("%.2f", Z)), size = 2.5, family = "Arial") +
  scale_fill_gradient2(low = "#1F78A8", mid = "white", high = "#C23B4A", midpoint = 0) +
  labs(title = "Top EPAS1 KO response genes across cohorts", x = NULL, y = NULL, fill = "Network Z") +
  theme_pub + theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(angle = 25, hjust = 1))

# C. Expression–strength relationship is a technical sensitivity check.
expr_plot <- ggplot(all_dr, aes(x = mean_log1p_cpm, y = abs(Z), colour = direction, shape = significant)) +
  geom_point(alpha = .8, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = .45, colour = "#16324F") +
  facet_wrap(~dataset, nrow = 1) + scale_colour_manual(values = pal) +
  labs(title = "Baseline expression versus perturbation strength",
       x = "Mean endothelial log1p(CP10K)", y = "Absolute network Z", colour = NULL, shape = "FDR < 0.05") + theme_pub

combined <- (volc / rank_plot / expr_plot) + plot_annotation(
  title = "Supplementary gene-level audit of EPAS1 virtual knockout",
  subtitle = "Gene-level displays support reproducibility assessment; program-level summaries remain the primary mechanistic evidence.",
  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 15),
                plot.subtitle = element_text(family = "Arial", size = 9, colour = "#4B5563")))
ggsave(file.path(out, "Figure_S26_EPAS1_KO_gene_level_audit.pdf"), combined, width = 12, height = 13, device = cairo_pdf)
png(file.path(out, "Figure_S26_EPAS1_KO_gene_level_audit.png"), width = 3600, height = 3900, res = 300)
print(combined)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
