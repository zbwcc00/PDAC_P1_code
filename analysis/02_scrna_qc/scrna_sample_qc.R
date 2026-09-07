options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

root <- "T:/"
base <- file.path(root, "data/02_scrna")
out_dir <- file.path(root, "analysis/02_scrna_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_files <- list.files(base, pattern = "matrix\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (!length(matrix_files)) stop("No matrix.mtx.gz files found")

read_lines <- function(path) {
  fread(path, header = FALSE, data.table = FALSE, showProgress = FALSE)
}

qc <- lapply(matrix_files, function(mtx_path) {
  parent <- dirname(mtx_path)
  prefix <- sub("_matrix\\.mtx\\.gz$", "", basename(mtx_path))
  barcode_path <- file.path(parent, paste0(prefix, "_barcodes.tsv.gz"))
  feature_path <- file.path(parent, paste0(prefix, "_features.tsv.gz"))
  if (!file.exists(feature_path)) feature_path <- file.path(parent, paste0(prefix, "_genes.tsv.gz"))
  if (!file.exists(barcode_path)) barcode_path <- file.path(parent, "barcodes.tsv.gz")
  if (!file.exists(feature_path)) feature_path <- file.path(parent, "features.tsv.gz")
  if (!file.exists(feature_path)) feature_path <- file.path(parent, "genes.tsv.gz")
  if (!file.exists(barcode_path) || !file.exists(feature_path)) {
    return(data.frame(dataset = NA, sample_label = NA, matrix_path = mtx_path, status = "missing_barcode_or_feature", cells = NA, genes = NA, median_umi = NA, median_genes_detected = NA))
  }
  matrix_obj <- readMM(gzfile(mtx_path))
  barcodes <- read_lines(barcode_path)
  features <- read_lines(feature_path)
  sample_label <- sub("_matrix\\.mtx\\.gz$", "", basename(mtx_path))
  if (grepl("GSE155698", mtx_path, fixed = TRUE)) {
    sample_label <- sub("\\.tar$", "", basename(dirname(dirname(mtx_path))))
  }
  dataset <- sub(".*(GSE[0-9]+).*", "\\1", mtx_path)
  data.frame(
    dataset = dataset,
    sample_label = sample_label,
    matrix_path = mtx_path,
    status = "ok",
    cells = ncol(matrix_obj),
    genes = nrow(matrix_obj),
    barcode_rows = nrow(barcodes),
    feature_rows = nrow(features),
    median_umi = as.numeric(median(Matrix::colSums(matrix_obj))),
    median_genes_detected = as.numeric(median(Matrix::colSums(matrix_obj > 0))),
    stringsAsFactors = FALSE
  )
})

qc <- rbindlist(qc, fill = TRUE)
fwrite(qc, file.path(out_dir, "PDAC_scRNA_sample_QC.tsv"), sep = "\t")
summary <- qc[, .(samples = .N, total_cells = sum(cells, na.rm = TRUE), min_cells = min(cells, na.rm = TRUE), median_cells = median(cells, na.rm = TRUE), max_cells = max(cells, na.rm = TRUE)), by = dataset]
fwrite(summary, file.path(out_dir, "PDAC_scRNA_dataset_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
