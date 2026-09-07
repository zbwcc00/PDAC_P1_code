options(stringsAsFactors = FALSE)
root <- "D:/PDAC_P1"
out <- file.path(root, "analysis/15_priority_figures/figure2_supplement")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2); library(patchwork)
  library(pheatmap); library(ggrepel); library(clusterProfiler); library(org.Hs.eg.db); library(ragg); library(svglite)
})
theme_pub <- theme_classic(base_family = "Arial", base_size = 8) +
  theme(axis.line = element_line(linewidth = .35, colour = "black"),
        axis.ticks = element_line(linewidth = .35, colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = 9, colour = "#111111"),
        plot.subtitle = element_text(size = 7, colour = "#444444"),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6),
        panel.grid = element_blank())
navy <- "#1F3B5B"; teal <- "#138A8A"; coral <- "#C85C5C"
save_pub <- function(p, stem, width = 12, height = 8) {
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "Arial"); print(p); dev.off()
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 300, type = "cairo-png"); print(p); dev.off()
  grDevices::tiff(paste0(stem, ".tiff"), width = width, height = height, units = "in", res = 300, compression = "lzw", type = "cairo"); print(p); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = width, height = height); print(p); dev.off()
}
read_geo_series <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)[1]
  titles <- strsplit(sub("^!Sample_title\\t", "", lines[grep("^!Sample_title", lines)[1]]), "\\t")[[1]]
  titles <- gsub('^"|"$', "", titles)
  accessions <- strsplit(sub("^!Sample_geo_accession\\t", "", lines[grep("^!Sample_geo_accession", lines)[1]]), "\\t")[[1]]
  accessions <- gsub('^"|"$', "", accessions)
  dat <- fread(path, skip = begin - 1, data.table = TRUE); setnames(dat, 1, "gene")
  expr <- as.matrix(dat[, -1, with = FALSE]); rownames(expr) <- dat$gene
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]; colnames(expr) <- accessions
  meta <- data.frame(accession = accessions, title = titles)
  meta$class <- ifelse(grepl("Primary-Pancreas", titles), "primary_pancreas",
    ifelse(grepl("Normal-Pancreas", titles), "normal_pancreas",
      ifelse(grepl("Normal-Vessel", titles), "normal_vessel",
        ifelse(grepl("Met-", titles), "metastasis", "other"))))
  list(expr = expr, meta = meta)
}
run_de <- function(expr, meta, paired = FALSE) {
  if (paired) {
    meta <- meta[meta$group %in% c("tumor", "adjacent") & !is.na(meta$pair_id), ]
    meta$group <- factor(meta$group, levels = c("adjacent", "tumor")); meta$pair_id <- factor(meta$pair_id)
    design <- model.matrix(~ pair_id + group, meta); mat <- expr[, meta$sample, drop = FALSE]; coef <- "grouptumor"
  } else {
    meta <- meta[meta$class %in% c("primary_pancreas", "normal_pancreas"), ]; meta$condition <- factor(meta$class, levels = c("normal_pancreas", "primary_pancreas"))
    design <- model.matrix(~ condition, meta); mat <- expr[, meta$accession, drop = FALSE]; coef <- "conditionprimary_pancreas"
  }
  fit <- eBayes(lmFit(mat, design)); tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); rownames(tt) <- NULL; list(tt = tt, meta = meta, expr = mat)
}
candidate <- fread(file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
candidate_genes <- unique(candidate$gene_id)

g624 <- fread(file.path(root, "analysis/08_program_validation/GSE62452/GSE62452_RMA_gene_symbol.tsv"))
g624_expr <- as.matrix(g624[, -1, with = FALSE]); rownames(g624_expr) <- g624[[1]]
g624_meta <- as.data.frame(fread(file.path(root, "analysis/08_program_validation/GSE62452/GSE62452_sample_metadata.tsv")))
de624 <- run_de(g624_expr, g624_meta, paired = TRUE)
g717 <- read_geo_series(file.path(root, "data/01_bulk/GSE71729/raw/GSE71729_series_matrix.txt.gz"))
de717 <- run_de(g717$expr, g717$meta, paired = FALSE)

pca_panel <- function(expr, meta, group_col, title, palette) {
  vars <- apply(expr, 1, var, na.rm = TRUE); keep <- names(sort(vars, decreasing = TRUE))[seq_len(min(1000, length(vars)))]
  pc <- prcomp(t(expr[keep, , drop = FALSE]), scale. = TRUE); v <- pc$sdev^2 / sum(pc$sdev^2)
  d <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], sample = rownames(pc$x), group = meta[[group_col]])
  ggplot(d, aes(PC1, PC2, colour = group)) + geom_point(size = 1.8, alpha = .85) +
    scale_colour_manual(values = palette, na.value = "#999999") +
    labs(title = title, x = sprintf("PC1 (%.1f%%)", 100*v[1]), y = sprintf("PC2 (%.1f%%)", 100*v[2]), colour = NULL) + theme_pub
}
volcano <- function(tt, title) {
  d <- tt; d$neglog10 <- -log10(pmax(d$P.Value, .Machine$double.xmin))
  d$flag <- ifelse(d$adj.P.Val < .05 & abs(d$logFC) >= 1, "FDR < 0.05 and |logFC| ≥ 1", "Other")
  lab <- d[d$gene %in% candidate_genes & d$adj.P.Val < .05, ]
  ggplot(d, aes(logFC, neglog10, colour = flag)) + geom_point(size = .65, alpha = .55) +
    scale_colour_manual(values = c("Other"="#B9C0C7", "FDR < 0.05 and |logFC| ≥ 1"=coral)) +
    geom_text_repel(data = lab, aes(label = gene), size = 2.4, max.overlaps = 12, colour = navy, min.segment.length = 0) +
    labs(title = title, x = "log2 fold change", y = "−log10(P value)", colour = NULL) + theme_pub
}
p624 <- pca_panel(de624$expr, de624$meta, "group", "GSE62452: paired bulk PCA", c(adjacent="#A7A7A7", tumor=navy))
p717 <- pca_panel(de717$expr, de717$meta, "class", "GSE71729: primary versus normal PCA", c(normal_pancreas="#A7A7A7", primary_pancreas=teal))
pvol624 <- volcano(de624$tt, "GSE62452 paired differential expression")
pvol717 <- volcano(de717$tt, "GSE71729 primary versus normal differential expression")

heat_genes <- intersect(candidate_genes, rownames(de717$expr)); heat_genes <- head(heat_genes, 30)
hm <- t(scale(t(de717$expr[heat_genes, de717$meta$accession, drop = FALSE]))); ord <- order(de717$meta$class)
ann_col <- data.frame(Group = de717$meta$class[ord]); rownames(ann_col) <- de717$meta$accession[ord]; hm <- hm[, ord, drop = FALSE]
png(file.path(out, "Figure_S24E_GSE71729_candidate_heatmap.png"), width = 3000, height = 2400, res = 300)
pheatmap(hm, annotation_col = ann_col, show_colnames = FALSE, fontsize_row = 7,
  color = colorRampPalette(c("#3B4F7C", "white", "#C85C5C"))(101), main = "Stable candidate genes in GSE71729")
dev.off()

make_gsea <- function(tt, title) {
  ids <- bitr(tt$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  m <- merge(tt[, c("gene", "t")], ids, by.x = "gene", by.y = "SYMBOL"); m <- m[!duplicated(m$ENTREZID), ]
  ranks <- m$t; names(ranks) <- m$ENTREZID; ranks <- sort(ranks, decreasing = TRUE)
  g <- tryCatch(gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE), error = function(e) NULL)
  if (is.null(g) || nrow(as.data.frame(g)) == 0) return(ggplot() + theme_void() + ggtitle(paste(title, "(no enriched terms)")))
  d <- head(as.data.frame(g)[order(as.data.frame(g)$p.adjust, -abs(as.data.frame(g)$NES)), ], 12); d$Description <- factor(d$Description, levels = rev(d$Description))
  ggplot(d, aes(NES, Description, size = -log10(p.adjust + 1e-12), colour = NES)) + geom_point() +
    scale_colour_gradient2(low = "#3B4F7C", mid = "white", high = coral, midpoint = 0) +
    labs(title = title, x = "Normalized enrichment score", y = NULL, size = "−log10(FDR)") + theme_pub
}
pgsea624 <- make_gsea(de624$tt, "GSE62452 GO biological-process GSEA")
pgsea717 <- make_gsea(de717$tt, "GSE71729 GO biological-process GSEA")
bulk_fig <- (p624 | p717) / (pvol624 | pvol717) / (pgsea624 | pgsea717) + plot_annotation(tag_levels = "A")
save_pub(bulk_fig, file.path(out, "Figure_S24_bulk_QC_DE_GSEA"), width = 13.5, height = 15)

sp <- fread(file.path(root, "analysis/09_spatial_validation/GSE282302/spatial_ecology/GSE282302_spatial_ecology_spots.tsv"))
score_cols <- c("malignant", "CAF", "myeloid", "T_NK", "endothelial"); sp <- sp[complete.cases(sp[, ..score_cols])]
set.seed(20260905); plot_sp <- sp[sample.int(.N, min(.N, 60000))]
pcsp <- prcomp(plot_sp[, ..score_cols], scale. = TRUE)
pcd <- data.frame(PC1 = pcsp$x[,1], PC2 = pcsp$x[,2], EPAS1_group = plot_sp$EPAS1_group)
psp_pca <- ggplot(pcd, aes(PC1, PC2, colour = EPAS1_group)) + geom_point(size = .25, alpha = .35) +
  scale_colour_manual(values = c(low="#3B4F7C", middle="#A7A7A7", high=coral), na.value = "#999999") +
  labs(title = "GSE282302 spatial ecology PCA", subtitle = "60,000 spots sampled for rendering; domain assignment used all eligible spots", colour = "EPAS1 group") + theme_pub
km <- kmeans(sp[, ..score_cols], centers = 4, nstart = 25); sp$domain <- factor(km$cluster)
rep_sample <- sp[, .N, by = sample][order(-N)]$sample[1]; rep <- sp[sample == rep_sample]
psp_domain <- ggplot(rep, aes(pxl_col_in_fullres, pxl_row_in_fullres, colour = domain)) + geom_point(size = .32, alpha = .7) +
  scale_y_reverse() + scale_colour_manual(values = c("#3B4F7C", "#138A8A", "#C58A24", "#C85C5C")) + coord_fixed() +
  labs(title = paste("GSE282302 spatial domains:", rep_sample), x = "Image x", y = "Image y", colour = "Domain") + theme_pub
sp297 <- fread(file.path(root, "analysis/09_spatial_validation/GSE297144/GSE297144_spot_anchor_scores.tsv"))
sp297 <- sp297[is.finite(EPAS1_logCPM) & is.finite(vascular_marker_score)]
set.seed(20260905); plot_sp297 <- sp297[sample.int(.N, min(.N, 60000))]
psp_bivar <- ggplot(plot_sp297, aes(EPAS1_logCPM, vascular_marker_score)) + geom_hex(bins = 45) + scale_fill_gradient(low = "#EEF2F6", high = navy) +
  geom_smooth(method = "lm", colour = coral, linewidth = .55, se = FALSE) +
  labs(title = "EPAS1–vascular-marker co-localization", subtitle = "GSE297144 spot-level logCPM; vascular score excludes EPAS1", x = "EPAS1 logCPM", y = "Vascular-marker score", fill = "Spot count") + theme_pub
pat282 <- fread(file.path(root, "analysis/09_spatial_validation/GSE282302/GSE282302_patient_spatial_anchor_results.tsv"))
pat297 <- fread(file.path(root, "analysis/09_spatial_validation/GSE297144/GSE297144_patient_spatial_anchor_results.tsv"))
sdat <- rbindlist(list(pat282[, .(dataset="GSE282302", patient, rho=median_EPAS1_vascular_rho)], pat297[, .(dataset="GSE297144", patient, rho=median_EPAS1_vascular_rho)]), fill = TRUE)
sdat$dataset <- factor(sdat$dataset, levels = c("GSE282302", "GSE297144"))
psp_patient <- ggplot(sdat, aes(rho, reorder(patient, rho), colour = dataset)) + geom_vline(xintercept = 0, linetype = 2, colour = "#777777") + geom_point(size = 2.2) +
  facet_wrap(~dataset, scales = "free_y") + scale_colour_manual(values = c("GSE282302"=navy, "GSE297144"=teal)) +
  labs(title = "Patient-level spatial replication", x = "Median Spearman ρ: EPAS1 versus vascular markers", y = NULL) + theme_pub + theme(legend.position = "none")
spatial_fig <- (psp_pca | psp_domain) / (psp_bivar | psp_patient) + plot_annotation(tag_levels = "A")
save_pub(spatial_fig, file.path(out, "Figure_S25_spatial_domains_colocalization"), width = 13.5, height = 10)

audit_lines <- c(
  "# Figure 2 supplemental data audit", "",
  paste0("- GSE62452: ", ncol(g624_expr), " arrays; ", nrow(g624_expr), " gene rows; paired limma used tumor–adjacent pairs."),
  paste0("- GSE71729: ", ncol(g717$expr), " samples; primary pancreas = ", sum(g717$meta$class == "primary_pancreas"), "; normal pancreas = ", sum(g717$meta$class == "normal_pancreas"), "."),
  paste0("- GSE282302: ", nrow(sp), " eligible spots across ", length(unique(sp$patient)), " patients; k-means domains fitted on all eligible spots."),
  "- Spatial PCA rendering samples 60,000 spots for plotting; all eligible spots were used for domain fitting and patient-level summaries.",
  "- GSE297144 contributes eight patient-level spatial anchor estimates; spots are not treated as independent patient replicates.",
  "- RCTD/cell2location is not promoted because a validated reference and comparable abundance output are not yet available; ecology scores are labeled as such."
)
writeLines(audit_lines, file.path(out, "Figure2_supplement_data_audit.md")); writeLines(capture.output(sessionInfo()), file.path(out, "40_R_sessionInfo.txt"))
message("Figure 2 supplemental outputs written to ", out)
