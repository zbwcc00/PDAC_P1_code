options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1"
outs <- file.path(root, "analysis/11_methodology_review")
audit_one <- function(gse, fn) {
  lines <- readLines(gzfile(fn), n = 2000)
  hdr <- lines[grep("^!Sample_geo_accession", lines)[1]]
  acc <- gsub('^"|"$', "", strsplit(sub("^!Sample_geo_accession\\t", "", hdr), "\\t")[[1]])
  title_line <- lines[grep("^!Sample_title", lines)[1]]
  titles <- gsub('^"|"$', "", strsplit(sub("^!Sample_title\\t", "", title_line), "\\t")[[1]])
  chars <- lines[grep("^!Sample_characteristics_ch1", lines)]
  data.frame(gse=gse, n_samples=length(acc), n_characteristics=length(chars), sample_accessions=paste(head(acc,5),collapse=","), titles_head=paste(head(titles,5),collapse=" || "), characteristic_fields=paste(sub("^!Sample_characteristics_ch1\\t", "", chars), collapse=" || "), stringsAsFactors=FALSE)
}
res <- rbind(audit_one("GSE21501", file.path(root,"data/01_bulk/GSE21501/raw/GSE21501_series_matrix_full.txt.gz")), audit_one("GSE62165", file.path(root,"data/01_bulk/GSE62165/raw/GSE62165_series_matrix_full.txt.gz")))
fwrite(res, file.path(outs,"prognostic_cohort_audit.tsv"), sep="\t")
sink(file.path(outs,"prognostic_cohort_audit.md")); cat("# Prognostic cohort audit\\n\\n"); print(res); sink()
