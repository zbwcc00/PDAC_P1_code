options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(SingleCellExperiment); library(slingshot); library(data.table)})
root <- "D:/PDAC_P1"
in_dir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")
out_dir <- file.path(root, "analysis/10_communication_pseudotime/Pseudotime")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(in_dir, pattern = "primary_tumor\\.rds$", full.names = TRUE)
candidate <- fread(file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
endo_up <- candidate[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id]
endo_down <- candidate[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id]
markers <- c("KDR","ESAM","PECAM1","VWF","EMCN","RAMP2","PLVAP","ENG","CA4","EPAS1","CTHRC1")
zrow <- function(x) as.numeric(scale(x))
signed <- function(m, up, down) {
  up <- intersect(up, rownames(m)); down <- intersect(down, rownames(m)); v <- rep(0, ncol(m)); n <- 0
  if (length(up)) {v <- v + rowMeans(apply(m[up,,drop=FALSE], 1, zrow)); n <- n+1}
  if (length(down)) {v <- v - rowMeans(apply(m[down,,drop=FALSE], 1, zrow)); n <- n+1}
  if (!n) rep(NA_real_, ncol(m)) else v/n
}
for (f in files) {
  dataset <- sub("__primary_tumor\\.rds$", "", basename(f))
  out_file <- file.path(out_dir, paste0(dataset, "_endothelial_pseudotime.tsv"))
  if (file.exists(out_file)) next
  obj <- readRDS(f); md <- obj@meta.data
  cells <- rownames(md)[md$major_lineage == "endothelial" & !is.na(md$patient_id)]
  if (length(cells) < 30) {fwrite(data.table(dataset=dataset, status="insufficient_cells", n_cells=length(cells)), out_file, sep="\t"); next}
  set.seed(20260902)
  if (length(cells) > 1200) cells <- sample(cells, 1200)
  mat <- GetAssayData(obj, assay="RNA", layer="data")[, cells, drop=FALSE]
  vars <- apply(as.matrix(mat), 1, var)
  gene_pool <- names(sort(vars, decreasing=TRUE))[seq_len(min(500, sum(vars > 0)))]
  gene_pool <- intersect(gene_pool, rownames(mat))
  x <- t(as.matrix(mat[gene_pool, , drop=FALSE]))
  pca <- prcomp(x, center=TRUE, scale.=TRUE, rank.=min(15, ncol(x), nrow(x)-1))
  pcs <- pca$x[, seq_len(min(10, ncol(pca$x))), drop=FALSE]
  k <- max(2, min(4, floor(nrow(pcs)/30)))
  km <- kmeans(pcs, centers=k, nstart=30)$cluster
  score <- signed(mat, endo_up, endo_down)
  cl_means <- tapply(score, km, mean, na.rm=TRUE)
  root <- names(which.min(cl_means))
  sce <- SingleCellExperiment(assays=list(logcounts=t(x)))
  reducedDims(sce)$PCA <- pcs
  colData(sce)$cluster <- factor(km)
  sds <- tryCatch(slingshot(sce, clusterLabels="cluster", start.clus=root, stretch=0), error=function(e) NULL)
  if (is.null(sds)) {pt <- as.numeric(scale(pcs[,1])); status <- "fallback_PC1"} else {pt <- slingPseudotime(sds)[,1]; status <- "slingshot"}
  out <- data.table(dataset=dataset, cell_id=rownames(pcs), patient_id=md[rownames(pcs), "patient_id"],
                    cluster=as.character(km), endothelial_program=score, pseudotime=pt, method=status)
  fwrite(out, out_file, sep="\t")
  cor_rows <- lapply(intersect(c(endo_up, endo_down, markers), rownames(mat)), function(g) data.table(gene=g, spearman_rho=cor(as.numeric(mat[g, cells]), pt, method="spearman", use="complete.obs")))
  fwrite(rbindlist(cor_rows), file.path(out_dir, paste0(dataset, "_pseudotime_gene_correlations.tsv")), sep="\t")
}
all_files <- list.files(out_dir, pattern="_endothelial_pseudotime.tsv$", full.names=TRUE)
all <- rbindlist(lapply(all_files, fread), fill=TRUE)
fwrite(all, file.path(out_dir, "all_endothelial_pseudotime.tsv"), sep="\t")
sink(file.path(out_dir, "endothelial_pseudotime_report.md")); cat("# Endothelial pseudotime\\n\\n"); cat("Trajectories were built independently per primary-tumor dataset using PCA on endothelial cells and Slingshot; root was chosen as the cluster with the lowest signed endothelial program score. This represents a transcriptional continuum, not proven lineage differentiation.\\n\\n"); print(all[, .(n=.N, n_patients=uniqueN(patient_id), program_rho=cor(endothelial_program,pseudotime,use="complete.obs")), by=dataset]); sink()
