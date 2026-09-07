options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(data.table)})
root <- "D:/PDAC_P1"
in_dir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")
out_dir <- file.path(root, "analysis/10_communication_pseudotime/Target_screen")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cand <- fread(file.path(root, "analysis/08_program_validation/endothelial_target_prioritization.tsv"))
genes <- unique(cand$gene_id[cand$target_priority %in% c("Tier1", "Tier2")])
files <- list.files(in_dir, pattern="primary_tumor\\.rds$", full.names=TRUE)
rows <- list()
for (f in files) {
  dataset <- sub("__primary_tumor\\.rds$", "", basename(f)); obj <- readRDS(f); md <- obj@meta.data
  avail <- intersect(genes, rownames(obj)); dat <- GetAssayData(obj, assay="RNA", layer="data")[avail,,drop=FALSE]
  for (ct in c("endothelial", "malignant_epithelial", "fibroblast_CAF", "myeloid", "T_NK")) {
    cells <- rownames(md)[md$major_lineage == ct]
    if (!length(cells)) next
    subm <- dat[, cells, drop=FALSE]
    for (g in avail) {
      vals <- as.numeric(subm[g,]); by_patient <- tapply(vals, md[cells,"patient_id"], mean, na.rm=TRUE)
      rows[[length(rows)+1]] <- data.table(dataset=dataset, cell_type=ct, gene=g, n_cells=length(cells), detection_fraction=mean(vals>0, na.rm=TRUE), median_expression=median(vals,na.rm=TRUE), mean_expression=mean(vals,na.rm=TRUE), n_patients_expressed=sum(by_patient>0,na.rm=TRUE), n_patients=length(unique(md[cells,"patient_id"])))
    }
  }
}
cov <- rbindlist(rows)
fwrite(cov, file.path(out_dir, "candidate_target_cell_coverage.tsv"), sep="\t")
endo <- cov[cell_type=="endothelial", .(datasets=uniqueN(dataset), mean_detection=mean(detection_fraction), mean_patients_expressed=mean(n_patients_expressed), total_cells=sum(n_cells)), by=gene]
spec <- dcast(cov, gene ~ cell_type, value.var="detection_fraction", fun.aggregate=mean)
fwrite(endo[order(-mean_detection)], file.path(out_dir, "endothelial_target_coverage_summary.tsv"), sep="\t")
fwrite(spec, file.path(out_dir, "candidate_target_detection_specificity.tsv"), sep="\t")
sink(file.path(out_dir, "target_cell_coverage_report.md")); cat("# Candidate target cell coverage\\n\\n"); cat("Detection fraction is calculated on normalized single-cell expression; patient-level expression support is reported to avoid cell-count pseudo-replication.\\n\\n"); print(endo[order(-mean_detection)]); cat("\\nDetection specificity matrix:\\n"); print(spec); sink()
