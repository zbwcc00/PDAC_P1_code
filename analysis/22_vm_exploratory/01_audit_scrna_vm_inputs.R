datasets <- c(
  GSE154778 = "D:/第二篇大论文/analysis/03_scrna_annotation/seurat_annotated/GSE154778__primary_tumor.rds",
  GSE155698 = "D:/第二篇大论文/analysis/03_scrna_annotation/seurat_annotated/GSE155698__primary_tumor.rds",
  GSE212966 = "D:/第二篇大论文/analysis/03_scrna_annotation/seurat_annotated/GSE212966__primary_tumor.rds"
)

out <- "C:/Users/HUAWEI/Documents/ChatGPT/Try/vm_exploratory_output"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

rows <- list()
genes <- c("EPAS1", "CDH5", "EPHA2", "MCAM", "ENG", "KDR", "FLT1", "PECAM1",
           "VWF", "EMCN", "LAMC2", "LAMA4", "COL4A1", "COL4A2", "MMP2", "MMP9",
           "ITGB1", "RHOA", "ROCK1", "ROCK2", "FAP", "COL1A1", "KRT19", "EPCAM")

for (dataset in names(datasets)) {
  object <- readRDS(datasets[[dataset]])
  metadata <- object[[]]
  gene_names <- rownames(object)
  for (field in colnames(metadata)) {
    values <- metadata[[field]]
    if (is.factor(values) || is.character(values)) {
      examples <- paste(head(sort(unique(as.character(values))), 20), collapse = " | ")
    } else {
      examples <- paste(head(sort(unique(values)), 10), collapse = " | ")
    }
    rows[[length(rows) + 1]] <- data.frame(
      dataset = dataset,
      n_cells = ncol(object),
      metadata_field = field,
      n_unique = length(unique(values)),
      examples = examples,
      stringsAsFactors = FALSE
    )
  }
  present <- data.frame(
    dataset = dataset,
    gene = genes,
    present = genes %in% gene_names,
    stringsAsFactors = FALSE
  )
  write.table(present, file.path(out, paste0(dataset, "_vm_gene_availability.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

write.table(do.call(rbind, rows), file.path(out, "scrna_input_metadata_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
