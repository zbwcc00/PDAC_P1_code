options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
})

root <- "T:/"
out_dir <- file.path(root, "analysis/01_tcga_qc")
counts <- readRDS(file.path(out_dir, "TCGA-PAAD_primary_tumor_counts_unstranded.rds"))
sample_qc <- fread(file.path(out_dir, "TCGA-PAAD_primary_tumor_qc.tsv"), data.table = FALSE)

keep <- rowSums(cpm(counts) > 1) >= 10
filtered_counts <- counts[keep, , drop = FALSE]
design <- matrix(1, nrow = ncol(filtered_counts), ncol = 1)
colnames(design) <- "intercept"
rownames(design) <- colnames(filtered_counts)
 y <- DGEList(counts = filtered_counts)
y <- calcNormFactors(y, method = "TMM")
log_cpm <- cpm(y, log = TRUE, prior.count = 2)

saveRDS(y, file.path(out_dir, "TCGA-PAAD_filtered_DGEList_TMM.rds"), compress = "xz")
saveRDS(log_cpm, file.path(out_dir, "TCGA-PAAD_filtered_logCPM.rds"), compress = "xz")
fwrite(data.table(gene_id = rownames(filtered_counts)), file.path(out_dir, "TCGA-PAAD_filtered_gene_ids.tsv"), sep = "\t")

pca <- prcomp(t(log_cpm), center = TRUE, scale. = FALSE)
pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
fwrite(as.data.table(pca_df), file.path(out_dir, "TCGA-PAAD_PCA_coordinates.tsv"), sep = "\t")

png(file.path(out_dir, "TCGA-PAAD_library_size_QC.png"), width = 1800, height = 1200, res = 180)
print(ggplot(sample_qc, aes(x = reorder(patient_id, library_size), y = library_size / 1e6)) + geom_col(fill = "#2C7FB8") + coord_flip() + labs(x = NULL, y = "Library size (million reads)", title = "TCGA-PAAD primary tumor library sizes") + theme_bw(base_size = 10))
dev.off()

png(file.path(out_dir, "TCGA-PAAD_PCA_QC.png"), width = 1800, height = 1400, res = 180)
print(ggplot(pca_df, aes(PC1, PC2)) + geom_point(size = 2.5, color = "#D95F02") + geom_text(aes(label = sample_id), size = 2.2, check_overlap = TRUE, vjust = -0.7) + labs(title = "TCGA-PAAD PCA of filtered TMM logCPM", x = sprintf("PC1 (%.1f%%)", pca_var[1]), y = sprintf("PC2 (%.1f%%)", pca_var[2])) + theme_bw(base_size = 11))
dev.off()

qc_summary <- data.frame(metric = c("input_genes", "filtered_genes", "patients", "filter_rule"), value = c(nrow(counts), nrow(filtered_counts), ncol(filtered_counts), "CPM > 1 in at least 10 patients"))
fwrite(as.data.table(qc_summary), file.path(out_dir, "TCGA-PAAD_normalization_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "normalization_sessionInfo.txt"))
message("Filtered genes: ", nrow(filtered_counts), "; patients: ", ncol(filtered_counts))
