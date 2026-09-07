options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
})

source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- pdac_paths()$project_root
manifest_path <- file.path(root, "analysis/02_scrna_qc/PDAC_scRNA_sample_manifest.tsv")
out_dir <- file.path(root, "analysis/02_scrna_qc")
checkpoint_dir <- file.path(out_dir, "seurat_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(manifest_path)
manifest <- manifest[matrix_status == "readable"]
resolve_manifest_path <- function(path, project_root) {
  normalized <- gsub("\\\\", "/", path)
  normalized <- sub("^[A-Za-z]:/第二篇大论文", project_root, normalized)
  normalized <- sub("^T:", project_root, normalized)
  normalized
}
manifest[, matrix_path := resolve_manifest_path(matrix_path, root)]
manifest[, barcode_path := resolve_manifest_path(barcode_path, root)]
manifest[, feature_path := resolve_manifest_path(feature_path, root)]

read_feature_names <- function(path) {
  x <- fread(path, header = FALSE, data.table = FALSE, showProgress = FALSE)
  if (ncol(x) >= 2) ids <- as.character(x[[2]]) else ids <- as.character(x[[1]])
  ids[is.na(ids) | ids == ""] <- as.character(x[[1]])[is.na(ids) | ids == ""]
  make.unique(ids)
}
read_barcodes <- function(path) {
  as.character(fread(path, header = FALSE, data.table = FALSE, showProgress = FALSE)[[1]])
}
read_matrix <- function(path) {
  Matrix::readMM(gzfile(path))
}

thresholds_for <- function(nfeature, ncount, pct_mt) {
  med_g <- median(nfeature); mad_g <- mad(nfeature, constant = 1.4826)
  med_u <- median(ncount); mad_u <- mad(ncount, constant = 1.4826)
  lower_gene <- max(100, floor(med_g - 3 * mad_g))
  lower_umi <- max(200, floor(med_u - 3 * mad_u))
  mt_q99 <- as.numeric(quantile(pct_mt, 0.99, names = FALSE, na.rm = TRUE))
  mt_med <- median(pct_mt, na.rm = TRUE)
  mt_mad <- mad(pct_mt, constant = 1.4826, na.rm = TRUE)
  mt_cut <- min(25, max(10, mt_med + 3 * mt_mad, mt_q99))
  upper_gene <- as.numeric(quantile(nfeature, 0.995, names = FALSE))
  upper_umi <- as.numeric(quantile(ncount, 0.995, names = FALSE))
  list(lower_gene = lower_gene, lower_umi = lower_umi, mt_cut = mt_cut,
       upper_gene_review = upper_gene, upper_umi_review = upper_umi,
       method = "lower genes/UMI=max(floor(median-3*MAD),minimum); mt=min(25,max(10,median+3*MAD,q99)); upper tails review-only")
}

qc_rows <- list(); object_by_context <- list(); errors <- list()
for (i in seq_len(nrow(manifest))) {
  rec <- manifest[i]
  message(sprintf("[%d/%d] %s %s", i, nrow(manifest), rec$dataset, rec$gsm))
  obj <- tryCatch({
    mat <- read_matrix(rec$matrix_path)
    genes <- read_feature_names(rec$feature_path)
    barcodes <- read_barcodes(rec$barcode_path)
    if (nrow(mat) != length(genes) || ncol(mat) != length(barcodes)) stop("dimension mismatch")
    rownames(mat) <- genes; colnames(mat) <- make.unique(barcodes)
    mat <- as(mat, "dgCMatrix")
    ncount <- Matrix::colSums(mat)
    nfeature <- Matrix::colSums(mat > 0)
    mt_idx <- grepl("^(MT-|MT\\.)", rownames(mat), ignore.case = TRUE)
    pct_mt <- if (any(mt_idx)) 100 * Matrix::colSums(mat[mt_idx, , drop = FALSE]) / pmax(ncount, 1) else rep(0, ncol(mat))
    th <- thresholds_for(nfeature, ncount, pct_mt)
    passes_core <- nfeature >= th$lower_gene & ncount >= th$lower_umi & pct_mt <= th$mt_cut

    sce <- SingleCellExperiment(assays = list(counts = mat))
    colData(sce)$sample_id <- rec$gsm
    dbl <- tryCatch({
      sce2 <- scDblFinder(sce, samples = rep(rec$gsm, ncol(sce)), verbose = FALSE)
      data.frame(doublet_score = as.numeric(colData(sce2)$scDblFinder.score),
                 doublet_class = as.character(colData(sce2)$scDblFinder.class),
                 doublet_status = "ok", stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(doublet_score = NA_real_, doublet_class = NA_character_, doublet_status = paste0("error:", conditionMessage(e)), stringsAsFactors = FALSE)
    })
    meta <- data.frame(cell = colnames(mat), nCount_RNA = as.numeric(ncount), nFeature_RNA = as.numeric(nfeature),
                       percent_mt = as.numeric(pct_mt), passes_core_QC = passes_core,
                       doublet_score = dbl$doublet_score, doublet_class = dbl$doublet_class,
                       doublet_status = dbl$doublet_status, dataset = rec$dataset, gsm = rec$gsm,
                       tissue_context = rec$tissue_context, patient_id = rec$patient_id,
                       stringsAsFactors = FALSE)
    qc_rows <<- append(qc_rows, list(data.frame(
      dataset = rec$dataset, gsm = rec$gsm, title = rec$title, tissue_context = rec$tissue_context,
      patient_id = rec$patient_id, cells_input = ncol(mat), genes_input = nrow(mat),
      median_umi = median(ncount), median_genes = median(nfeature), median_percent_mt = median(pct_mt),
      lower_gene_cutoff = th$lower_gene, lower_umi_cutoff = th$lower_umi, mt_cutoff = th$mt_cut,
      upper_gene_review = th$upper_gene_review, upper_umi_review = th$upper_umi_review,
      cells_pass_core = sum(passes_core), cells_flagged_core = sum(!passes_core),
      doublets_called = sum(dbl$doublet_class == "doublet", na.rm = TRUE), doublet_rate = mean(dbl$doublet_class == "doublet", na.rm = TRUE),
      doublet_status = paste(unique(dbl$doublet_status), collapse = ";"), threshold_method = th$method,
      stringsAsFactors = FALSE)))
    obj <- CreateSeuratObject(counts = mat, project = rec$dataset, min.cells = 0, min.features = 0)
    rownames(meta) <- meta$cell
    obj <- AddMetaData(obj, meta = meta[, setdiff(names(meta), "cell"), drop = FALSE])
    obj
  }, error = function(e) {
    errors <<- append(errors, list(data.frame(dataset = rec$dataset, gsm = rec$gsm, error = conditionMessage(e), stringsAsFactors = FALSE)))
    NULL
  })
  if (!is.null(obj)) {
    key <- paste(rec$dataset, rec$tissue_context, sep = "__")
    if (is.null(object_by_context[[key]])) object_by_context[[key]] <- list()
    object_by_context[[key]][[rec$gsm]] <- obj
  }
  rm(obj); gc(verbose = FALSE)
}

qc <- rbindlist(qc_rows, fill = TRUE)
fwrite(qc, file.path(out_dir, "PDAC_scRNA_cell_QC_summary.tsv"), sep = "\t", na = "NA")
if (length(errors)) fwrite(rbindlist(errors, fill = TRUE), file.path(out_dir, "PDAC_scRNA_QC_errors.tsv"), sep = "\t", na = "NA")

context_summary <- qc[, .(samples = .N, cells_input = sum(cells_input), cells_pass_core = sum(cells_pass_core),
                          median_pass_fraction = median(cells_pass_core / cells_input), median_umi = median(median_umi),
                          median_genes = median(median_genes), median_mt = median(median_percent_mt),
                          median_doublet_rate = median(doublet_rate, na.rm = TRUE)), by = .(dataset, tissue_context)]
fwrite(context_summary, file.path(out_dir, "PDAC_scRNA_context_QC_summary.tsv"), sep = "\t", na = "NA")

for (key in names(object_by_context)) {
  pieces <- object_by_context[[key]]
  if (!length(pieces)) next
  merged <- if (length(pieces) == 1) pieces[[1]] else merge(pieces[[1]], y = pieces[-1], add.cell.ids = names(pieces), merge.data = FALSE)
  merged@misc$qc_policy <- list(core_metrics = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
                                thresholds = "sample-specific; see PDAC_scRNA_cell_QC_summary.tsv",
                                doublet_caller = "scDblFinder 1.20.2 per sample", raw_counts_preserved = TRUE)
  safe_key <- gsub("[^A-Za-z0-9_]+", "_", key)
  saveRDS(merged, file.path(checkpoint_dir, paste0(safe_key, ".rds")), compress = "xz")
  rm(merged); gc(verbose = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "scrna_qc_seurat_sessionInfo.txt"))
print(context_summary)
