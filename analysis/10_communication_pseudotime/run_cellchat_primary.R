options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(CellChat); library(data.table)})
source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- pdac_paths()$project_root
in_dir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")
out_dir <- file.path(root, "analysis/10_communication_pseudotime/CellChat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = 8 * 1024^3)
future::plan("sequential")

files <- list.files(in_dir, pattern = "primary_tumor\\.rds$", full.names = TRUE)
for (f in files) {
  dataset <- sub("__primary_tumor\\.rds$", "", basename(f))
  out_rds <- file.path(out_dir, paste0(dataset, "_cellchat.rds"))
  if (file.exists(out_rds)) next
  message("Processing ", dataset)
  obj <- readRDS(f)
  md <- obj@meta.data
  keep_types <- c("malignant_epithelial", "epithelial", "fibroblast_CAF", "endothelial", "myeloid", "T_NK", "CD8", "B_cell", "mast")
  keep_cells <- rownames(md)[md$major_lineage %in% keep_types & !is.na(md$patient_id)]
  md2 <- md[keep_cells, , drop = FALSE]
  set.seed(20260902)
  sampled <- unlist(lapply(split(rownames(md2), interaction(md2$patient_id, md2$major_lineage, drop = TRUE)), function(x) sample(x, min(length(x), 80))))
  sampled <- intersect(sampled, colnames(obj))
  md2 <- md[sampled, , drop = FALSE]
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")[, sampled, drop = FALSE]
  meta_chat <- data.frame(group = factor(md2$major_lineage), row.names = rownames(md2))
  cellchat <- createCellChat(object = expr, meta = meta_chat, group.by = "group")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, out_rds, compress = "xz")
  comm <- subsetCommunication(cellchat)
  fwrite(as.data.table(comm), file.path(out_dir, paste0(dataset, "_all_communications.tsv")), sep = "\t")
  endo_in <- subsetCommunication(cellchat, targets.use = "endothelial")
  endo_out <- subsetCommunication(cellchat, sources.use = "endothelial")
  fwrite(as.data.table(endo_in), file.path(out_dir, paste0(dataset, "_to_endothelial.tsv")), sep = "\t")
  fwrite(as.data.table(endo_out), file.path(out_dir, paste0(dataset, "_from_endothelial.tsv")), sep = "\t")
}
sink(file.path(out_dir, "CellChat_run_report.md")); cat("# CellChat primary-tumor analysis\\n\\n"); cat("Each dataset was analyzed independently; up to 80 cells per patient × major lineage were sampled. Results are exploratory and require patient-level validation.\\n"); sink()
