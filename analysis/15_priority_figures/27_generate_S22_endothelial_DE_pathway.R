base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(clusterProfiler); library(org.Hs.eg.db)
})
theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "#555555"),
        axis.title = element_text(size = 9), axis.text = element_text(colour = "#222222"),
        legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "#F2F4F6", colour = NA),
        strip.text = element_text(face = "bold"))
save_pub <- function(p, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  png(file.path(out_dir, paste0(stem, ".png")), width = round(width * 300), height = round(height * 300), res = 300, type = "cairo-png")
  print(p); dev.off()
}

de_files <- c(
  GSE155698_edgeR = file.path(base_dir, "analysis/06_pseudobulk_DE/edgeR/GSE155698/endothelial_DE.tsv"),
  GSE155698_DESeq2 = file.path(base_dir, "analysis/06_pseudobulk_DE/DESeq2/GSE155698/endothelial_DE.tsv"),
  GSE212966_edgeR = file.path(base_dir, "analysis/06_pseudobulk_DE/edgeR/GSE212966/endothelial_DE.tsv"),
  GSE212966_DESeq2 = file.path(base_dir, "analysis/06_pseudobulk_DE/DESeq2/GSE212966/endothelial_DE.tsv")
)
read_de <- function(path, label) {
  x <- read.delim(path, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  if ("logFC" %in% names(x)) { x$effect <- x$logFC; x$fdr <- x$FDR }
  else { x$effect <- x$log2FoldChange; x$fdr <- x$padj }
  x$gene_id <- gsub('^"|"$', "", x$gene_id)
  x$test <- label
  x[, c("gene_id", "effect", "fdr", "test")]
}
de <- do.call(rbind, Map(read_de, de_files, names(de_files))); rownames(de) <- NULL
de$dataset <- ifelse(grepl("GSE155698", de$test), "GSE155698", "GSE212966")
de$method <- ifelse(grepl("edgeR", de$test), "edgeR", "DESeq2")

candidate_genes <- c("EPAS1", "TACC1", "MARCKS", "HERPUD1")
candidate <- de[de$gene_id %in% candidate_genes, ]
candidate$gene_id <- factor(candidate$gene_id, levels = rev(candidate_genes))
candidate$test <- factor(candidate$test, levels = names(de_files))
p22a <- ggplot(candidate, aes(test, gene_id, fill = effect)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", effect)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Endothelial candidate effects", x = NULL, y = NULL, fill = "log2 fold-change") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))

# Select genes with the most reproducible direction across all four endothelial tests.
wide <- reshape(de[, c("gene_id", "test", "effect")], idvar = "gene_id", timevar = "test", direction = "wide")
test_cols <- paste0("effect.", names(de_files))
test_cols <- test_cols[test_cols %in% names(wide)]
wide$direction_consistent <- apply(wide[, test_cols, drop = FALSE], 1, function(z) all(z > 0) || all(z < 0))
wide$mean_abs_effect <- rowMeans(abs(wide[, test_cols, drop = FALSE]), na.rm = TRUE)
shared <- wide[wide$direction_consistent & rowSums(!is.na(wide[, test_cols, drop = FALSE])) == length(test_cols), ]
shared <- shared[order(-shared$mean_abs_effect), , drop = FALSE]
top_genes <- unique(c(candidate_genes, head(shared$gene_id, 12)))
top_genes <- top_genes[top_genes %in% unique(de$gene_id)]
heat <- de[de$gene_id %in% top_genes, c("gene_id", "test", "effect")]
gene_order <- top_genes[order(match(top_genes, candidate_genes), top_genes)]
heat$gene_id <- factor(heat$gene_id, levels = rev(unique(gene_order)))
heat$test <- factor(heat$test, levels = names(de_files))
p22b <- ggplot(heat, aes(test, gene_id, fill = effect)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.1f", effect)), size = 2.25) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Cross-method DE consistency", x = NULL, y = NULL, fill = "log2 fold-change") +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7), axis.text.y = element_text(size = 7))

# Pathway gate: reuse the validated local enrichment table when available so
# figure regeneration is deterministic and does not depend on KEGG network access.
existing_all <- file.path(out_dir, "S22_endothelial_DE_GO_KEGG_all.tsv")
existing_recurrent <- file.path(out_dir, "S22_endothelial_DE_GO_KEGG_recurrent.tsv")
if (file.exists(existing_all) && file.exists(existing_recurrent)) {
  path_df <- read.delim(existing_all, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  path_sum <- read.delim(existing_recurrent, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
} else {
path_rows <- list(); k <- 0
for (test_name in names(de_files)) {
  x <- de[de$test == test_name & !is.na(de$fdr) & de$fdr < 0.05, ]
  for (direction in c("up", "down")) {
    genes <- if (direction == "up") x$gene_id[x$effect > 0] else x$gene_id[x$effect < 0]
    genes <- unique(genes[genes != "" & !is.na(genes)])
    if (length(genes) < 5) next
    ego <- tryCatch(enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", readable = FALSE), error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      z <- as.data.frame(ego); z <- z[z$p.adjust < 0.05, , drop = FALSE]
      if (nrow(z) > 0) { z <- z[seq_len(min(30, nrow(z))), ]; k <- k + 1; path_rows[[k]] <- data.frame(source = "GO-BP", direction = direction, test = test_name, term = z$Description, padj = z$p.adjust, stringsAsFactors = FALSE) }
    }
    ids <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, error = function(e) character())
    if (length(unique(ids)) >= 5) {
      ek <- tryCatch(enrichKEGG(gene = unique(ids), organism = "hsa", pAdjustMethod = "BH"), error = function(e) NULL)
      if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
        z <- as.data.frame(ek); z <- z[z$p.adjust < 0.05, , drop = FALSE]
        if (nrow(z) > 0) { z <- z[seq_len(min(30, nrow(z))), ]; k <- k + 1; path_rows[[k]] <- data.frame(source = "KEGG", direction = direction, test = test_name, term = z$Description, padj = z$p.adjust, stringsAsFactors = FALSE) }
      }
    }
  }
}
path_df <- if (length(path_rows)) do.call(rbind, path_rows) else data.frame(source = character(), direction = character(), test = character(), term = character(), padj = numeric())
if (nrow(path_df)) {
  path_count <- aggregate(padj ~ source + direction + term, path_df, length)
  names(path_count)[names(path_count) == "padj"] <- "n_tests"
  path_min <- aggregate(padj ~ source + direction + term, path_df, min)
  names(path_min)[names(path_min) == "padj"] <- "min_padj"
  path_sum <- merge(path_count, path_min, by = c("source", "direction", "term"))
  path_sum <- path_sum[path_sum$n_tests >= 2, , drop = FALSE]
  path_sum <- path_sum[order(-path_sum$n_tests, path_sum$min_padj), , drop = FALSE]
  path_sum <- head(path_sum, 10)
} else path_sum <- data.frame()
}
write.table(path_df, file.path(out_dir, "S22_endothelial_DE_GO_KEGG_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(path_sum, file.path(out_dir, "S22_endothelial_DE_GO_KEGG_recurrent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

if (nrow(path_sum)) {
  path_sum$term <- factor(path_sum$term, levels = rev(unique(path_sum$term)))
  p22c <- ggplot(path_sum, aes(n_tests, term, colour = source, size = -log10(min_padj))) + geom_point() + scale_x_continuous(breaks = 2:4, limits = c(1.8, 4.2)) + scale_colour_manual(values = c("GO-BP" = "#009E73", "KEGG" = "#D55E00")) + labs(title = "Recurrent pathway enrichment", x = "Number of four DE tests", y = NULL, size = "−log10 minimum FDR") + theme_pub
} else {
  p22c <- ggplot() + annotate("text", x = 0, y = 0, label = "No GO-BP/KEGG term recurred\n(FDR < 0.05 in >=2 of 4 tests)", size = 4, colour = "#555555") + xlim(-1, 1) + ylim(-1, 1) + labs(title = "Recurrent pathway enrichment", x = NULL, y = NULL) + theme_void() + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
}
p22 <- (p22a | p22b) / p22c + plot_annotation(title = "Supplementary Figure S22. Endothelial DE and pathway-convergence gate", subtitle = "Tumour-versus-adjacent-normal endothelial pseudobulk; four tests (2 cohorts x 2 methods). Pathways retained only when FDR < 0.05 in at least 2 tests.") & theme(plot.title = element_text(face = "bold", size = 10.5), plot.subtitle = element_text(size = 7.2, margin = margin(t = 2, b = 3)))
save_pub(p22, "Figure_S22_endothelial_DE_pathway_gate", 9.2, 7.2)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "27_R_sessionInfo.txt"))
