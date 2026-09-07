options(stringsAsFactors = FALSE)

out_dir <- "D:/PDAC_P1/analysis/11_methodology_review"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_geo <- function(path, accession) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con)
  get_field <- function(prefix) {
    hit <- lines[startsWith(lines, prefix)]
    if (!length(hit)) return(character())
    vals <- strsplit(hit[1], "\t", fixed = TRUE)[[1]][-1]
    vals <- gsub('^"|"$', "", vals)
    vals
  }
  begin <- which(lines == "!series_matrix_table_begin")
  end <- which(lines == "!series_matrix_table_end")
  if (length(begin) != 1L || length(end) != 1L) stop("matrix table markers not found")
  tab <- read.delim(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                    check.names = FALSE, stringsAsFactors = FALSE,
                    quote = "\"", comment.char = "")
  sample_ids <- get_field("!Sample_geo_accession")
  sample_titles <- get_field("!Sample_title")
  characteristics <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  characteristic_keys <- vapply(strsplit(characteristics, "\t", fixed = TRUE), function(x) {
    x <- x[-1]
    x <- gsub('^"|"$', "", x)
    key <- sub(":.*$", "", x)
    paste(unique(key), collapse = "|")
  }, character(1))
  sample_n <- ncol(tab) - 1L
  meta <- data.frame(
    accession = accession,
    sample_id = if (length(sample_ids)) sample_ids else colnames(tab)[-1],
    title = if (length(sample_titles)) sample_titles else NA_character_,
    stringsAsFactors = FALSE
  )
  meta$characteristic_keys <- if (length(characteristic_keys)) paste(characteristic_keys, collapse = ";") else NA_character_
  meta$n_samples_matrix <- sample_n
  meta$n_features <- nrow(tab)
  meta$platform <- get_field("!Series_platform_id")[1] %||% NA_character_
  pheno <- data.frame(sample_id = meta$sample_id, stringsAsFactors = FALSE)
  if (length(characteristics)) {
    for (line in characteristics) {
      vals <- gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]][-1])
      keys <- sub(":.*$", "", vals)
      values <- sub("^[^:]*:\\s*", "", vals)
      for (k in unique(keys)) {
        if (!nzchar(k)) next
        col <- make.names(k)
        if (!col %in% names(pheno)) pheno[[col]] <- NA_character_
        pheno[[col]] <- ifelse(keys == k, values, pheno[[col]])
      }
    }
  }
  list(meta = meta, pheno = pheno, tab = tab, lines = lines)
}

`%||%` <- function(x, y) if (length(x) && !is.na(x[1]) && nzchar(x[1])) x[1] else y

datasets <- list(
  GSE21501 = "D:/PDAC_P1/data/01_bulk/GSE21501/raw/GSE21501_series_matrix_full.txt.gz",
  GSE62165 = "D:/PDAC_P1/data/01_bulk/GSE62165/raw/GSE62165_series_matrix_full.txt.gz"
)

all_meta <- list()
all_pheno <- list()
summary <- list()
for (acc in names(datasets)) {
  x <- parse_geo(datasets[[acc]], acc)
  all_meta[[acc]] <- x$meta
  x$pheno$accession <- acc
  all_pheno[[acc]] <- x$pheno
  lines <- x$lines
  series_title <- sub('^!Series_title\t', '', lines[startsWith(lines, "!Series_title")][1])
  series_title <- gsub('^"|"$', "", series_title)
  characteristics <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  characteristic_text <- paste(head(unique(unlist(lapply(strsplit(characteristics, "\t", fixed = TRUE), function(z) gsub('^"|"$', "", z[-1])))), 30), collapse = " | ")
  summary[[acc]] <- data.frame(
    accession = acc,
    title = series_title,
    n_samples = ncol(x$tab) - 1L,
    n_features = nrow(x$tab),
    feature_id_column = colnames(x$tab)[1],
    platform = x$meta$platform[1],
    sample_characteristic_examples = characteristic_text,
    stringsAsFactors = FALSE
  )
}

meta_out <- do.call(rbind, all_meta)
all_pheno_names <- unique(unlist(lapply(all_pheno, names)))
all_pheno <- lapply(all_pheno, function(z) {
  z[setdiff(all_pheno_names, names(z))] <- NA_character_
  z[, all_pheno_names, drop = FALSE]
})
pheno_out <- do.call(rbind, all_pheno)
summary_out <- do.call(rbind, summary)
write.table(meta_out, file.path(out_dir, "bulk_geo_sample_metadata_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pheno_out, file.path(out_dir, "bulk_geo_phenotype_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(summary_out, file.path(out_dir, "bulk_geo_matrix_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

field_hits <- names(pheno_out)[grepl("surv|os|dfs|pfs|death|censor|follow|month|day|year", names(pheno_out), ignore.case = TRUE)]
field_hits <- setdiff(field_hits, c("accession", "sample_id"))

report <- c(
  "# GEO bulk matrix audit",
  "",
  "This audit is descriptive and does not treat array intensities as RNA-seq counts.",
  "",
  paste0("- GSE21501: ", summary_out$n_samples[summary_out$accession == "GSE21501"], " samples; ", summary_out$n_features[summary_out$accession == "GSE21501"], " probe rows."),
  paste0("- GSE62165: ", summary_out$n_samples[summary_out$accession == "GSE62165"], " samples; ", summary_out$n_features[summary_out$accession == "GSE62165"], " probe rows."),
  "- Both are expression-profiling-by-array series matrices; use limma after probe annotation and cohort-specific normalization.",
  "- Survival variables must be confirmed from GEO supplementary phenotype tables or the original publication; absence in the series matrix is recorded as missing, not imputed.",
  paste0("- Phenotype fields matching survival/time/status keywords: ", if (length(field_hits)) paste(field_hits, collapse = ", ") else "none detected."),
  "",
  "Outputs: bulk_geo_matrix_audit.tsv, bulk_geo_sample_metadata_audit.tsv, and bulk_geo_phenotype_table.tsv"
)
writeLines(report, file.path(out_dir, "bulk_geo_audit_report.md"))
