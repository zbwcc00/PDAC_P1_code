options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(limma)})
root <- "D:/PDAC_P1"
in_file <- file.path(root, "data/01_bulk/GSE71729/raw/GSE71729_series_matrix.txt.gz")
out_dir <- file.path(root, "analysis/08_program_validation/GSE71729")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(in_file)) stop("Missing series matrix: ", in_file)

lines <- readLines(gzfile(in_file))
title_line <- lines[grep("^!Sample_title", lines)[1]]
titles <- sub("^!Sample_title\\t", "", title_line)
titles <- strsplit(titles, "\\t")[[1]]
titles <- gsub('^"|"$', "", titles)
acc_line <- lines[grep("^!Sample_geo_accession", lines)[1]]
accessions <- strsplit(sub("^!Sample_geo_accession\\t", "", acc_line), "\\t")[[1]]
accessions <- gsub('^"|"$', "", accessions)
dat <- fread(in_file, skip = grep("^!series_matrix_table_begin", lines)[1] - 1, data.table = TRUE)
setnames(dat, 1, "gene")
expr <- as.matrix(dat[, -1, with = FALSE]); rownames(expr) <- dat$gene
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
colnames(expr) <- accessions
meta <- data.frame(accession = accessions, title = titles, stringsAsFactors = FALSE)
meta$group_detail <- sub("^[0-9]+-", "", meta$title)
meta$sample_class <- ifelse(grepl("Primary-Pancreas", meta$title), "primary_pancreas",
                      ifelse(grepl("Normal-Pancreas", meta$title), "normal_pancreas",
                      ifelse(grepl("Normal-Vessel", meta$title), "normal_vessel",
                      ifelse(grepl("Met-", meta$title), "metastasis", "other"))))
fwrite(meta, file.path(out_dir, "GSE71729_sample_metadata.tsv"), sep = "\t")

candidate <- fread(file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
endo_up <- candidate[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id]
endo_down <- candidate[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id]
tnk_up <- candidate[major_lineage == "T_NK" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id]
tnk_down <- candidate[major_lineage == "T_NK" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id]
zrow <- function(x) as.numeric(scale(x))
signed <- function(m, up, down) {
  up <- intersect(up, rownames(m)); down <- intersect(down, rownames(m)); v <- rep(0, ncol(m)); n <- 0
  if (length(up)) {v <- v + rowMeans(apply(m[up,,drop=FALSE], 1, zrow)); n <- n+1}
  if (length(down)) {v <- v - rowMeans(apply(m[down,,drop=FALSE], 1, zrow)); n <- n+1}
  if (!n) rep(NA_real_, ncol(m)) else v/n
}
scores <- data.frame(accession = colnames(expr), endothelial = signed(expr, endo_up, endo_down),
                     T_NK = signed(expr, tnk_up, tnk_down), stringsAsFactors = FALSE)
scores <- merge(scores, meta, by = "accession", sort = FALSE)
fwrite(scores, file.path(out_dir, "GSE71729_program_scores.tsv"), sep = "\t")

run_contrast <- function(df, g1, g0, label) {
  d <- df[df$sample_class %in% c(g1, g0), ]; d$condition <- factor(d$sample_class, levels = c(g0, g1))
  if (nrow(d) < 6 || length(unique(d$condition)) < 2) return(data.frame())
  design <- model.matrix(~ condition, d)
  coef_name <- colnames(design)[2]
  out <- lapply(c("endothelial", "T_NK"), function(p) {
    fit <- eBayes(lmFit(matrix(d[[p]], nrow = 1), design))
    tt <- topTable(fit, coef = coef_name, number = 1, sort.by = "none")
    data.frame(contrast = label, program = p, n_g1 = sum(d$condition == g1), n_g0 = sum(d$condition == g0),
               logFC = tt$logFC, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
  })
  rbindlist(out)
}
res <- rbindlist(list(run_contrast(scores, "primary_pancreas", "normal_pancreas", "primary_vs_normal"),
                      run_contrast(scores, "normal_vessel", "normal_pancreas", "vessel_vs_normal"),
                      run_contrast(scores, "metastasis", "normal_pancreas", "met_vs_normal")), fill = TRUE)
fwrite(res, file.path(out_dir, "GSE71729_program_contrasts.tsv"), sep = "\t")
primary <- scores[scores$sample_class %in% c("primary_pancreas", "normal_pancreas"), ]
primary$condition <- factor(primary$sample_class, levels = c("normal_pancreas", "primary_pancreas"))
design_primary <- model.matrix(~ condition, primary)
candidate_genes <- unique(c(endo_up, endo_down, tnk_up, tnk_down))
gene_keep <- intersect(candidate_genes, rownames(expr))
gene_fit <- eBayes(lmFit(expr[gene_keep, primary$accession, drop = FALSE], design_primary))
gene_tt <- topTable(gene_fit, coef = colnames(design_primary)[2], number = Inf, sort.by = "none")
gene_tt$gene <- rownames(gene_tt); rownames(gene_tt) <- NULL
gene_tt <- gene_tt[, c("gene", "logFC", "P.Value", "adj.P.Val")]
fwrite(gene_tt, file.path(out_dir, "GSE71729_primary_vs_normal_candidate_gene_results.tsv"), sep = "\t")
summary <- data.table(scores)[, .(n = .N, endothelial_mean = mean(endothelial, na.rm=TRUE), T_NK_mean = mean(T_NK, na.rm=TRUE)), by = sample_class]
fwrite(summary, file.path(out_dir, "GSE71729_program_group_summary.tsv"), sep = "\t")
sink(file.path(out_dir, "GSE71729_validation_report.md")); cat("# GSE71729 external bulk validation\\n\\n"); cat("- Processed GEO series matrix (GPL20769), 357 samples; values are GEO-provided normalized expression.\\n"); cat("- Candidate program scores use gene-wise z-scoring and signed tumor-up minus tumor-down scoring.\\n\\n"); print(summary); cat("\\n"); print(res); sink()
