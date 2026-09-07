options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(oligo)
  library(limma)
  library(data.table)
  library(hugene10sttranscriptcluster.db)
  library(AnnotationDbi)
})

root <- "D:/PDAC_P1"
in_dir <- file.path(root, "data/01_bulk/GSE62452/raw/CEL")
out_dir <- file.path(root, "analysis/08_program_validation/GSE62452")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cel <- sort(list.files(in_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE))
if (length(cel) < 100) stop("Expected 100+ CEL files; found ", length(cel))

eset_file <- file.path(out_dir, "GSE62452_RMA_core_eset.rds")
if (file.exists(eset_file)) {
  message("Loading cached RMA object")
  eset <- readRDS(eset_file)
} else {
  message("Reading ", length(cel), " CEL files")
  raw <- read.celfiles(cel, verbose = TRUE)
  eset <- rma(raw, target = "core", normalize = TRUE)
  saveRDS(eset, eset_file, compress = "xz")
}
expr_probe <- exprs(eset)
colnames(expr_probe) <- sub("\\.CEL\\.gz$", "", basename(colnames(expr_probe)), ignore.case = TRUE)

probe_ids <- rownames(expr_probe)
ann <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                             keys = probe_ids,
                             keytype = "PROBEID",
                             columns = c("SYMBOL"))
ann <- ann[!is.na(ann$SYMBOL) & ann$SYMBOL != "", c("PROBEID", "SYMBOL")]
ann <- ann[!duplicated(ann$PROBEID), ]
keep <- intersect(probe_ids, ann$PROBEID)
expr_probe <- expr_probe[keep, , drop = FALSE]
ann <- ann[match(rownames(expr_probe), ann$PROBEID), ]
expr_gene <- limma::avereps(expr_probe, ID = ann$SYMBOL)
write.table(as.data.frame(expr_gene), file.path(out_dir, "GSE62452_RMA_gene_symbol.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

meta <- data.frame(sample = colnames(expr_gene), stringsAsFactors = FALSE)
meta$gsm <- sub("_.*$", "", meta$sample)
meta$group <- ifelse(grepl("_T$", meta$sample), "tumor",
                     ifelse(grepl("_N$", meta$sample), "adjacent", "other"))
meta$title <- sub("^[^_]+_[0-9]+_", "", meta$sample)
meta$pair_id <- ifelse(grepl("^E[0-9]+", meta$title), sub("^(E[0-9]+).*", "\\1", meta$title), NA)
meta$pair_id[meta$group != "tumor" & meta$group != "adjacent"] <- NA
fwrite(meta, file.path(out_dir, "GSE62452_sample_metadata.tsv"), sep = "\t")

candidate <- fread(file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
endo <- candidate[major_lineage == "endothelial" & stable_high_confidence == TRUE, gene_id]
tnk <- candidate[major_lineage == "T_NK" & stable_high_confidence == TRUE, gene_id]
endo_up <- candidate[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id]
endo_down <- candidate[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id]
tnk_up <- candidate[major_lineage == "T_NK" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id]
tnk_down <- candidate[major_lineage == "T_NK" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id]

zscore <- function(x) as.numeric(scale(x))
score_signed <- function(mat, up, down) {
  up <- intersect(up, rownames(mat)); down <- intersect(down, rownames(mat))
  if (!length(up) && !length(down)) return(rep(NA_real_, ncol(mat)))
  val <- rep(0, ncol(mat)); n <- 0
  if (length(up)) { val <- val + rowMeans(apply(mat[up, , drop = FALSE], 1, zscore)); n <- n + 1 }
  if (length(down)) { val <- val - rowMeans(apply(mat[down, , drop = FALSE], 1, zscore)); n <- n + 1 }
  val / n
}
scores <- data.frame(sample = colnames(expr_gene),
                     endothelial = score_signed(expr_gene, endo_up, endo_down),
                     T_NK = score_signed(expr_gene, tnk_up, tnk_down),
                     stringsAsFactors = FALSE)
scores <- merge(scores, meta, by = "sample", sort = FALSE)
fwrite(scores, file.path(out_dir, "GSE62452_program_scores.tsv"), sep = "\t")

paired <- scores[!is.na(scores$pair_id) & scores$group %in% c("tumor", "adjacent"), ]
paired <- paired[!duplicated(paired[, c("pair_id", "group")]), ]
paired <- paired[paired$pair_id %in% names(which(table(paired$pair_id) == 2)), ]
paired$group <- factor(paired$group, levels = c("adjacent", "tumor"))
paired$pair_id <- factor(paired$pair_id)
fit_program <- function(y) {
  fit <- lmFit(matrix(y, nrow = 1), model.matrix(~ pair_id + group, paired))
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "grouptumor", number = 1, sort.by = "none")
  data.frame(logFC = tt$logFC, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
}
prog_res <- rbind(endothelial = fit_program(paired$endothelial), T_NK = fit_program(paired$T_NK))
prog_res$program <- rownames(prog_res); rownames(prog_res) <- NULL
prog_res$n_pairs <- length(unique(paired$pair_id))
fwrite(prog_res, file.path(out_dir, "GSE62452_paired_program_results.tsv"), sep = "\t")

gene_rows <- list()
for (g in intersect(c(endo, tnk), rownames(expr_gene))) {
  fit <- lmFit(matrix(expr_gene[g, paired$sample], nrow = 1), model.matrix(~ pair_id + group, paired))
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "grouptumor", number = 1, sort.by = "none")
  gene_rows[[g]] <- data.frame(gene = g, logFC = tt$logFC, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
}
gene_res <- rbindlist(gene_rows, fill = TRUE)
fwrite(gene_res, file.path(out_dir, "GSE62452_paired_candidate_gene_results.tsv"), sep = "\t")

sink(file.path(out_dir, "GSE62452_validation_report.md"))
cat("# GSE62452 paired external validation\\n\\n")
cat("- Platform: GPL6244 Affymetrix Human Gene 1.0 ST; RMA core summarization and probe-to-symbol averaging.\\n")
cat("- Total CEL files:", ncol(expr_gene), "\\n")
cat("- Paired samples:", length(unique(paired$pair_id)), "\\n")
cat("- Endothelial genes available:", length(intersect(endo, rownames(expr_gene))), "/", length(endo), "\\n")
cat("- T/NK genes available:", length(intersect(tnk, rownames(expr_gene))), "/", length(tnk), "\\n\\n")
print(prog_res)
sink()
message("Done")
