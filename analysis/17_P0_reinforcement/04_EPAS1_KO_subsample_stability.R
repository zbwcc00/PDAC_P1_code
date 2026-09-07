options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(scTenifoldKnk); library(data.table)})

root <- "D:/PDAC_P1"
outdir <- file.path(root, "analysis/17_P0_reinforcement")
objdir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")
datasets <- c(
  GSE155698_primary = "GSE155698__primary_tumor.rds",
  GSE212966_primary = "GSE212966__primary_tumor.rds"
)
programs <- list(
  endothelial_downstream = c("CLDN5", "KDR", "PECAM1", "VWF", "EMCN", "ESAM", "RAMP2", "PLVAP", "ENG"),
  angiogenesis_downstream = c("VEGFA", "KDR", "PECAM1", "VWF", "EMCN", "ENG")
)
rows <- list()
for (dataset in names(datasets)) {
  reference <- readRDS(file.path(root, "analysis/13_virtual_perturbation/results", paste0(dataset, "__endothelial__EPAS1__KO.rds")))
  genes <- rownames(reference$tensorNetworks$WT)
  object <- readRDS(file.path(objdir, datasets[[dataset]]))
  assay <- if ("RNA" %in% names(object@assays)) "RNA" else DefaultAssay(object)
  counts <- GetAssayData(object, assay = assay, layer = "counts")
  metadata <- object@meta.data
  cells <- rownames(metadata)[metadata$major_lineage == "endothelial"]
  sub <- counts[intersect(genes, rownames(counts)), cells, drop = FALSE]
  sub <- sub[, Matrix::colSums(sub) >= 100, drop = FALSE]
  sub <- sub[genes[genes %in% rownames(sub)], , drop = FALSE]
  for (replicate_id in seq_len(5L)) {
    set.seed(20261000 + replicate_id)
    fit <- scTenifoldKnk(sub, qc = FALSE, gKO = "EPAS1", nc_nNet = 2,
      nc_nCells = min(100, ncol(sub)), nc_nComp = 3, td_K = 3, ma_nDim = 2, nCores = 1)
    dr <- fit$diffRegulation
    rows[[length(rows) + 1L]] <- rbindlist(lapply(names(programs), function(endpoint) {
      d <- dr[match(intersect(programs[[endpoint]], dr$gene), dr$gene), ]
      data.table(
        dataset = dataset, replicate_id = replicate_id, endpoint = endpoint,
        n_FDR05 = sum(d$p.adj < 0.05, na.rm = TRUE),
        downstream_network_score = sum(-log10(pmax(d$p.adj, 1e-300)), na.rm = TRUE),
        downstream_mean_abs_Z = mean(abs(d$Z), na.rm = TRUE)
      )
    }))
    rm(fit); gc()
  }
  rm(object, counts, sub, reference); gc()
}
results <- rbindlist(rows)
fwrite(results, file.path(outdir, "EPAS1_KO_positive_cohort_subsample_repeats.tsv"), sep = "\t")
summary <- results[, .(
  n_repeats = .N,
  median_network_score = median(downstream_network_score),
  min_network_score = min(downstream_network_score),
  max_network_score = max(downstream_network_score),
  positive_FDR_repeat_fraction = mean(n_FDR05 > 0)
), by = .(dataset, endpoint)]
fwrite(summary, file.path(outdir, "EPAS1_KO_positive_cohort_subsample_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(outdir, "04_sessionInfo.txt"))

