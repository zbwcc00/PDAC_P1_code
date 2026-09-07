options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

root <- "T:/"
registry <- file.path(root, "data/00_registry/GDC")
raw_dir <- file.path(root, "data/01_bulk/TCGA-PAAD/raw")
out_dir <- file.path(root, "analysis/01_tcga_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

map <- fread(file.path(registry, "TCGA-PAAD_primary_tumor_sample_map.tsv"), sep = "\t", data.table = FALSE)
sample_type_col <- "cases.0.samples.0.sample_type"
patient_col <- "cases.0.submitter_id"
sample_col <- "cases.0.samples.0.submitter_id"
map <- map[map[[sample_type_col]] == "Primary Tumor", , drop = FALSE]
map <- map[!duplicated(map[[patient_col]]), , drop = FALSE]
map$patient_id <- map[[patient_col]]
map$sample_id <- map[[sample_col]]
map$local_path <- file.path(raw_dir, map$file_id, map$file_name)
missing_files <- map[!file.exists(map$local_path), c("file_id", "file_name", patient_col), drop = FALSE]
fwrite(missing_files, file.path(out_dir, "missing_expression_files.tsv"), sep = "\t")
if (nrow(missing_files) > 0) stop("Missing local expression files: ", nrow(missing_files))

read_counts <- function(path) {
  tab <- fread(path, sep = "\t", skip = "gene_id", select = c("gene_id", "gene_name", "gene_type", "unstranded"), data.table = FALSE)
  tab <- tab[grepl("^ENSG", tab$gene_id), c("gene_id", "gene_name", "gene_type", "unstranded"), drop = FALSE]
  tab$gene_id <- sub("\\..*$", "", tab$gene_id)
  tab$unstranded <- as.numeric(tab$unstranded)
  tab <- tab[!is.na(tab$unstranded), , drop = FALSE]
  annotation <- tab[!duplicated(tab$gene_id), c("gene_id", "gene_name", "gene_type"), drop = FALSE]
  values <- aggregate(unstranded ~ gene_id, data = tab, FUN = sum)
  values$gene_name <- annotation$gene_name[match(values$gene_id, annotation$gene_id)]
  values$gene_type <- annotation$gene_type[match(values$gene_id, annotation$gene_id)]
  values[, c("gene_id", "gene_name", "gene_type", "unstranded"), drop = FALSE]
}

first <- read_counts(map$local_path[1])
gene_ids <- first$gene_id
gene_names <- first$gene_name[match(gene_ids, first$gene_id)]
gene_types <- first$gene_type[match(gene_ids, first$gene_id)]
counts <- matrix(0, nrow = length(gene_ids), ncol = nrow(map), dimnames = list(gene_ids, map[[patient_col]]))
read_qc <- vector("list", nrow(map))

for (i in seq_len(nrow(map))) {
  tab <- if (i == 1) first else read_counts(map$local_path[i])
  index <- match(gene_ids, tab$gene_id)
  counts[, i] <- ifelse(is.na(index), 0, tab$unstranded[index])
  read_qc[[i]] <- data.frame(patient_id = map[[patient_col]][i], sample_id = map[[sample_col]][i], file_id = map$file_id[i], file_name = map$file_name[i], genes_detected = sum(tab$unstranded > 0), library_size = sum(tab$unstranded, na.rm = TRUE), mapped_to_reference = sum(!is.na(index)))
  if (i %% 20 == 0 || i == nrow(map)) message("Read ", i, "/", nrow(map))
}

gene_annotation <- data.frame(gene_id = gene_ids, gene_name = gene_names, gene_type = gene_types)
sample_qc <- rbindlist(read_qc, fill = TRUE)
sample_metadata <- merge(map, sample_qc, by = c("patient_id", "sample_id", "file_id", "file_name"), all.x = TRUE, sort = FALSE)
saveRDS(counts, file.path(out_dir, "TCGA-PAAD_primary_tumor_counts_unstranded.rds"), compress = "xz")
saveRDS(gene_annotation, file.path(out_dir, "TCGA-PAAD_gene_annotation.rds"), compress = "xz")
fwrite(as.data.table(sample_metadata), file.path(out_dir, "TCGA-PAAD_primary_tumor_sample_metadata.tsv"), sep = "\t")
fwrite(as.data.table(sample_qc), file.path(out_dir, "TCGA-PAAD_primary_tumor_qc.tsv"), sep = "\t")
qc_summary <- data.frame(metric = c("patients", "genes", "min_library_size", "median_library_size", "max_library_size", "min_genes_detected", "median_genes_detected", "max_genes_detected"), value = c(ncol(counts), nrow(counts), min(sample_qc$library_size), median(sample_qc$library_size), max(sample_qc$library_size), min(sample_qc$genes_detected), median(sample_qc$genes_detected), max(sample_qc$genes_detected)))
fwrite(as.data.table(qc_summary), file.path(out_dir, "TCGA-PAAD_qc_summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Completed: ", ncol(counts), " patients x ", nrow(counts), " genes")
