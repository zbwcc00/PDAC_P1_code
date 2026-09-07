base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures/scrna_additional")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#555555"),
        axis.title = element_text(size = 9),
        axis.text = element_text(colour = "#222222"),
        legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "#F2F4F6", colour = NA),
        strip.text = element_text(face = "bold"))
save_pub <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  png(file.path(out_dir, paste0(stem, ".png")), width = round(width * 300), height = round(height * 300), res = 300, type = "cairo-png")
  print(p); dev.off()
}
read_tsv <- function(x) read.delim(x, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")

# 1. Sample-level QC and cellular composition.
qc <- read_tsv(file.path(base_dir, "analysis/02_scrna_qc/PDAC_scRNA_cell_QC_summary.tsv"))
qc$dataset <- gsub('"', '', qc$dataset)
qc$dataset <- factor(qc$dataset, levels = c("GSE154778", "GSE155698", "GSE212966"))
qc$pass_fraction <- as.numeric(qc$pass_fraction)
qc$median_genes <- as.numeric(qc$median_genes)
qc$median_umi <- as.numeric(qc$median_umi)
qc$median_percent_mt <- as.numeric(qc$median_percent_mt)
p_qc1 <- ggplot(qc, aes(dataset, pass_fraction, fill = dataset)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "#333333") +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.05)) +
  scale_fill_manual(values = c(GSE154778 = "#0072B2", GSE155698 = "#009E73", GSE212966 = "#D55E00")) +
  labs(title = "Core QC pass fraction", x = NULL, y = "Cells passing core QC") + theme_pub + theme(legend.position = "none")
p_qc2 <- ggplot(qc, aes(dataset, median_genes, fill = dataset)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "#333333") +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.7) +
  scale_fill_manual(values = c(GSE154778 = "#0072B2", GSE155698 = "#009E73", GSE212966 = "#D55E00")) +
  labs(title = "Median detected genes", x = NULL, y = "Genes per cell") + theme_pub + theme(legend.position = "none")
p_qc3 <- ggplot(qc, aes(dataset, median_umi, fill = dataset)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "#333333") +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.7) +
  scale_fill_manual(values = c(GSE154778 = "#0072B2", GSE155698 = "#009E73", GSE212966 = "#D55E00")) +
  labs(title = "Median UMI", x = NULL, y = "UMI per cell") + theme_pub + theme(legend.position = "none")
p_qc4 <- ggplot(qc, aes(dataset, median_percent_mt, fill = dataset)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "#333333") +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.7) +
  scale_fill_manual(values = c(GSE154778 = "#0072B2", GSE155698 = "#009E73", GSE212966 = "#D55E00")) +
  labs(title = "Mitochondrial fraction", x = NULL, y = "% mitochondrial reads") + theme_pub + theme(legend.position = "none")
save_pub((p_qc1 | p_qc2) / (p_qc3 | p_qc4) + plot_annotation(title = "Supplementary Figure S15. Single-cell sample-level quality control", subtitle = "Each point is one GEO sample; boxplots summarize the three discovery cohorts."), "Figure_S15_scRNA_sample_QC", 8.6, 6.6)

# 2. Marker DotPlot-like display aggregated by cohort and lineage.
dataset_files <- c(GSE154778 = "GSE154778__primary_tumor.rds", GSE155698 = "GSE155698__primary_tumor.rds", GSE212966 = "GSE212966__primary_tumor.rds")
marker_map <- list(
  endothelial = c("EMCN", "KDR", "ESAM", "PECAM1", "VWF"),
  fibroblast_CAF = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),
  myeloid = c("LST1", "TYROBP", "CTSS", "FCER1G", "SPP1"),
  T_NK = c("CD3D", "CD3E", "TRBC1", "NKG7", "CCL5"),
  B_cell = c("CD79A", "MS4A1", "CD74", "CD37"),
  epithelial = c("KRT8", "KRT18", "KRT19", "EPCAM"),
  acinar = c("PRSS1", "REG1A", "AMY2A"),
  mast = c("TPSAB1", "KIT")
)
markers <- unique(unlist(marker_map))
marker_df <- list()
for (ds in names(dataset_files)) {
  obj <- readRDS(file.path(base_dir, "analysis/03_scrna_annotation/seurat_annotated", dataset_files[[ds]]))
  avail <- intersect(markers, rownames(obj))
  dat <- FetchData(obj, vars = c(avail, "major_lineage"))
  for (gene in avail) {
    tmp <- data.frame(dataset = ds, gene = gene, major_lineage = dat$major_lineage, value = dat[[gene]])
    tmp <- tmp[!is.na(tmp$major_lineage), ]
    agg <- aggregate(value ~ gene + major_lineage, tmp, function(x) mean(x, na.rm = TRUE))
    pct <- aggregate(value ~ gene + major_lineage, tmp, function(x) mean(x > 0, na.rm = TRUE))
    names(agg)[3] <- "avg"; names(pct)[3] <- "pct"
    m <- merge(agg, pct, by = c("gene", "major_lineage")); m$dataset <- ds
    marker_df[[length(marker_df) + 1]] <- m
  }
}
marker_df <- do.call(rbind, marker_df)
marker_df <- marker_df[marker_df$major_lineage %in% names(marker_map), ]
marker_df$gene <- factor(marker_df$gene, levels = rev(markers))
marker_df$major_lineage <- factor(marker_df$major_lineage, levels = names(marker_map))
p_dot <- ggplot(marker_df, aes(major_lineage, gene, size = pct, colour = avg)) +
  geom_point(alpha = 0.9) + facet_wrap(~dataset, nrow = 1) +
  scale_size(range = c(0.5, 5), labels = scales::percent_format(accuracy = 25)) +
  scale_colour_gradient(low = "#F4F7FB", high = "#0072B2") +
  labs(title = "Canonical lineage markers across discovery cohorts", x = NULL, y = NULL, size = "Detected cells", colour = "Mean expression") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7), axis.text.y = element_text(size = 7), legend.position = "right")
save_pub(p_dot, "Figure_S16_scRNA_marker_DotPlot", 10.2, 7.0)

# 3. Endothelial pseudobulk differential expression: four prespecified tests.
de_files <- c(
  GSE155698_edgeR = file.path(base_dir, "analysis/06_pseudobulk_DE/edgeR/GSE155698/endothelial_DE.tsv"),
  GSE155698_DESeq2 = file.path(base_dir, "analysis/06_pseudobulk_DE/DESeq2/GSE155698/endothelial_DE.tsv"),
  GSE212966_edgeR = file.path(base_dir, "analysis/06_pseudobulk_DE/edgeR/GSE212966/endothelial_DE.tsv"),
  GSE212966_DESeq2 = file.path(base_dir, "analysis/06_pseudobulk_DE/DESeq2/GSE212966/endothelial_DE.tsv")
)
de_list <- lapply(names(de_files), function(nm) {
  d <- read_tsv(de_files[[nm]]); d$test <- nm
  d$logFC_use <- if ("logFC" %in% names(d)) as.numeric(d$logFC) else as.numeric(d$log2FoldChange)
  d$FDR_use <- if ("FDR" %in% names(d)) as.numeric(d$FDR) else as.numeric(d$padj)
  d[, c("gene_id", "logFC_use", "FDR_use", "test")]
})
de <- do.call(rbind, de_list)
de <- de[is.finite(de$logFC_use) & is.finite(de$FDR_use), ]
de$gene_id <- gsub('"', '', de$gene_id)
de$neglog10FDR <- -log10(pmax(de$FDR_use, 1e-300)); de$significant <- de$FDR_use < 0.05
label_genes <- c("EPAS1", "TACC1", "MARCKS", "HERPUD1", "PECAM1", "VWF", "KDR", "EMCN", "ESAM")
de$label <- ifelse(de$gene_id %in% label_genes, de$gene_id, NA)
p_vol <- ggplot(de, aes(logFC_use, neglog10FDR, colour = significant)) +
  geom_point(size = 0.65, alpha = 0.42) +
  geom_point(data = de[de$gene_id %in% label_genes, ], colour = "#D55E00", size = 1.4) +
  ggrepel::geom_text_repel(data = de[de$gene_id %in% label_genes, ], aes(label = gene_id), size = 2.8, max.overlaps = 20, colour = "#222222", box.padding = 0.3, seed = 20260905) +
  facet_wrap(~test, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = c(`FALSE` = "#B8C0C8", `TRUE` = "#0072B2"), labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  labs(title = "Endothelial pseudobulk differential expression", x = "Tumour versus reference log2 fold-change", y = "-log10 adjusted P", colour = "FDR < 0.05") + theme_pub + theme(legend.position = "bottom")
save_pub(p_vol, "Figure_S17_endothelial_pseudobulk_volcano", 9.2, 7.0)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "32_R_sessionInfo.txt"))
