options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
})

root <- "D:/PDAC_P1"
soft_dir <- file.path(root, "data/00_registry/GEO")
data_dir <- file.path(root, "data/02_scrna")
out_dir <- file.path(root, "analysis/02_scrna_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_soft <- function(path, dataset) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- which(grepl("^\\^SAMPLE", lines))
  if (!length(starts)) stop("No sample blocks in ", path)
  ends <- c(starts[-1] - 1L, length(lines))
  rows <- lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    accession <- sub("^\\^SAMPLE\\s*=\\s*", "", block[1])
    get_one <- function(tag) {
      hit <- block[grepl(paste0("^", tag, "\\s*=\\s*"), block)]
      if (!length(hit)) return(NA_character_)
      sub(paste0("^", tag, "\\s*=\\s*"), "", hit[1])
    }
    chars <- block[grepl("^!Sample_characteristics_ch1", block)]
    chars <- sub("^!Sample_characteristics_ch1\\s*=\\s*", "", chars)
    data.frame(
      dataset = dataset,
      gsm = accession,
      title = get_one("!Sample_title"),
      source_name = get_one("!Sample_source_name_ch1"),
      characteristics = paste(chars, collapse = " | "),
      platform = get_one("!Sample_platform_id"),
      stringsAsFactors = FALSE
    )
  })
  rbindlist(rows, fill = TRUE)
}

soft_files <- c(
  GSE154778 = file.path(soft_dir, "GSE154778_family.soft"),
  GSE155698 = file.path(soft_dir, "GSE155698_family.soft"),
  GSE212966 = file.path(soft_dir, "GSE212966_family.soft")
)
meta <- rbindlist(lapply(names(soft_files), function(ds) parse_soft(soft_files[[ds]], ds)), fill = TRUE)

meta[, tissue_context := "unknown"]
meta[dataset == "GSE154778" & grepl("Primary", title, ignore.case = TRUE), tissue_context := "primary_tumor"]
meta[dataset == "GSE154778" & grepl("Metastatic", title, ignore.case = TRUE), tissue_context := "metastasis"]
meta[dataset == "GSE155698" & grepl("PANCREAS TUMOR", source_name, ignore.case = TRUE), tissue_context := "primary_tumor"]
meta[dataset == "GSE155698" & grepl("Adjancent|Adjacent", source_name, ignore.case = TRUE), tissue_context := "adjacent_normal"]
meta[dataset == "GSE155698" & grepl("PBMC", source_name, ignore.case = TRUE) & grepl("Healthy", title, ignore.case = TRUE), tissue_context := "healthy_pbmc"]
meta[dataset == "GSE155698" & grepl("PBMC", source_name, ignore.case = TRUE) & !grepl("Healthy", title, ignore.case = TRUE), tissue_context := "pdac_pbmc"]
meta[dataset == "GSE212966" & grepl("^PDAC", title, ignore.case = TRUE), tissue_context := "primary_tumor"]
meta[dataset == "GSE212966" & grepl("^ADJ", title, ignore.case = TRUE), tissue_context := "adjacent_normal"]
meta[, treatment_status := ifelse(grepl("Treatment Na.?ve|No treatment", characteristics, ignore.case = TRUE), "untreated", "not_reported")]
meta[, patient_id := fifelse(dataset == "GSE154778", sub(":.*", "", title),
                             fifelse(dataset == "GSE155698", sub(".*_(\\d+[A-Z]?)$", "\\1", title),
                                     sub(",.*", "", title)))]
meta[dataset == "GSE155698", patient_id := sub("[A-Z]$", "", patient_id)]

matrix_files <- list.files(data_dir, pattern = "matrix\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
manifest <- rbindlist(lapply(matrix_files, function(path) {
  ds <- sub(".*(GSE[0-9]+).*", "\\1", path)
  acc <- sub("^(GSM[0-9]+).*", "\\1", basename(path))
  if (grepl("GSE155698", path, fixed = TRUE)) {
    acc <- sub(".*(GSM[0-9]+).*", "\\1", path)
  }
  parent <- dirname(path)
  barcode <- if (grepl("GSE155698", path, fixed = TRUE)) file.path(parent, "barcodes.tsv.gz") else sub("_matrix\\.mtx\\.gz$", "_barcodes.tsv.gz", path)
  feature <- if (grepl("GSE155698", path, fixed = TRUE)) file.path(parent, "features.tsv.gz") else sub("_matrix\\.mtx\\.gz$", "_features.tsv.gz", path)
  if (!file.exists(feature) && !grepl("GSE155698", path, fixed = TRUE)) feature <- sub("_matrix\\.mtx\\.gz$", "_genes.tsv.gz", path)
  if (!file.exists(feature) && grepl("GSE155698", path, fixed = TRUE)) feature <- file.path(parent, "genes.tsv.gz")
  data.frame(dataset = ds, gsm = acc, matrix_path = normalizePath(path, winslash = "/", mustWork = FALSE),
             barcode_path = normalizePath(barcode, winslash = "/", mustWork = FALSE),
             feature_path = normalizePath(feature, winslash = "/", mustWork = FALSE),
             matrix_exists = file.exists(path), barcode_exists = file.exists(barcode), feature_exists = file.exists(feature),
             stringsAsFactors = FALSE)
}), fill = TRUE)

out <- merge(meta, manifest, by = c("dataset", "gsm"), all = TRUE, sort = FALSE)
out[, matrix_status := fifelse(is.na(matrix_path), "not_found", fifelse(matrix_exists & barcode_exists & feature_exists, "readable", "incomplete_bundle"))]
setcolorder(out, c("dataset", "gsm", "title", "tissue_context", "patient_id", "source_name", "characteristics", "treatment_status", "platform", "matrix_status", "matrix_path", "barcode_path", "feature_path", "matrix_exists", "barcode_exists", "feature_exists"))
fwrite(out, file.path(out_dir, "PDAC_scRNA_sample_manifest.tsv"), sep = "\t", na = "NA")
fwrite(out, file.path(soft_dir, "PDAC_scRNA_sample_manifest.tsv"), sep = "\t", na = "NA")
summary <- out[, .(samples = .N, readable_bundles = sum(matrix_status == "readable", na.rm = TRUE),
                   primary_or_tumor = sum(tissue_context == "primary_tumor", na.rm = TRUE),
                   adjacent = sum(tissue_context == "adjacent_normal", na.rm = TRUE),
                   metastasis = sum(tissue_context == "metastasis", na.rm = TRUE),
                   pdac_pbmc = sum(tissue_context == "pdac_pbmc", na.rm = TRUE),
                   healthy_pbmc = sum(tissue_context == "healthy_pbmc", na.rm = TRUE),
                   unknown_context = sum(tissue_context == "unknown", na.rm = TRUE)), by = dataset]
fwrite(summary, file.path(out_dir, "PDAC_scRNA_manifest_summary.tsv"), sep = "\t", na = "NA")
print(summary)
