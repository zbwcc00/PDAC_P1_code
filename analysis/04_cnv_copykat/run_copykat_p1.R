options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table); library(copykat)})
set.seed(20260901)
root <- "D:/PDAC_P1"
qc_dir <- file.path(root, "analysis/02_scrna_qc/seurat_checkpoints")
out_dir <- file.path(root, "analysis/04_cnv_copykat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

configs <- list(
  GSE154778 = list(primary="GSE154778__primary_tumor.rds", reference="GSE154778__metastasis.rds", ref_context="metastasis"),
  GSE155698 = list(primary="GSE155698__primary_tumor.rds", reference="GSE155698__adjacent_normal.rds", ref_context="adjacent_normal"),
  GSE212966 = list(primary="GSE212966__primary_tumor.rds", reference="GSE212966__adjacent_normal.rds", ref_context="adjacent_normal")
)

extract_cells <- function(obj, max_per_sample=500) {
  if("RNA" %in% names(obj@assays)) obj <- JoinLayers(obj, assay="RNA")
  md <- as.data.frame(obj[[]]); md$cell <- rownames(md); setDT(md)
  md <- md[passes_core_QC %in% TRUE & (is.na(doublet_class) | doublet_class != "doublet")]
  if(!nrow(md)) stop("No QC-passing singlet cells")
  chosen <- md[, .SD[sample.int(.N, min(.N, max_per_sample))], by=.(gsm)]
  counts <- LayerData(obj, assay="RNA", layer="counts")[, chosen$cell, drop=FALSE]
  colnames(counts) <- chosen$cell
  list(counts=counts, metadata=chosen)
}

filter_genes <- function(mat) {
  keep <- Matrix::rowSums(mat > 0) >= 20
  keep <- keep & !grepl("^(MT-|RPL|RPS|MTRNR)", rownames(mat), ignore.case=TRUE)
  mat[keep, , drop=FALSE]
}

all_calls <- list(); run_meta <- list()
for(ds in names(configs)) {
  pred_name <- paste0("PDAC_", ds, "_copykat_prediction.txt")
  if(file.exists(file.path(out_dir, ds, pred_name))) {
    message("Skipping completed ", ds)
    next
  }
  cfg <- configs[[ds]]; message("Preparing ", ds)
  p_obj <- readRDS(file.path(qc_dir, cfg$primary)); r_obj <- readRDS(file.path(qc_dir, cfg$reference))
  p <- extract_cells(p_obj, 100); r <- extract_cells(r_obj, 50)
  p$metadata$cnv_role <- "tumor_or_metastasis"
  r$metadata$cnv_role <- "reference_context"
  common_genes <- intersect(rownames(p$counts), rownames(r$counts))
  if(length(common_genes) < 5000) stop("Too few common genes for CNV inference: ", length(common_genes))
  p$counts <- p$counts[common_genes, , drop=FALSE]
  r$counts <- r$counts[common_genes, , drop=FALSE]
  mat <- cbind(p$counts, r$counts)
  md <- rbind(p$metadata, r$metadata)
  rownames(md) <- md$cell
  mat <- filter_genes(mat)
  colnames(mat) <- make.unique(colnames(mat))
  rownames(md) <- colnames(mat)
  run_dir <- file.path(out_dir, ds); dir.create(run_dir, recursive=TRUE, showWarnings=FALSE)
  saveRDS(list(counts=mat, metadata=md), file.path(run_dir, "copykat_input_sampled.rds"), compress="xz")
  old <- getwd(); setwd(run_dir); on.exit(setwd(old), add=TRUE)
  message("Running CopyKAT ", ds, " cells=", ncol(mat), " genes=", nrow(mat))
  ck <- tryCatch(copykat(rawmat=as.matrix(mat), id.type="S", cell.line="no", ngene.chr=5,
                         min.gene.per.cell=200, LOW.DR=0.05, UP.DR=0.1, win.size=25,
                         sam.name=paste0("PDAC_", ds), output.seg="TRUE", plot.genes="FALSE",
                         genome="hg20", n.cores=1), error=function(e) e)
  if(inherits(ck,"error")) {
    fwrite(data.frame(dataset=ds,status="error",error=conditionMessage(ck)), file.path(run_dir,"copykat_status.tsv"), sep="\t")
    next
  }
  saveRDS(ck, file.path(run_dir, "copykat_result.rds"), compress="xz")
  pred <- ck$prediction
  if(!is.null(pred)) {
    pred <- as.data.table(pred); names(pred)[1:2] <- c("cell","copykat_prediction")
    pred <- merge(pred, md[,.(cell, dataset, gsm, tissue_context, patient_id, cnv_role)], by="cell", all.x=TRUE)
    fwrite(pred, file.path(run_dir,"copykat_cell_calls.tsv"), sep="\t", na="NA")
    all_calls[[ds]] <- pred
    run_meta[[ds]] <- data.frame(dataset=ds, status="ok", sampled_cells=ncol(mat), genes=nrow(mat),
      predicted_aneuploid=sum(grepl("aneuploid", pred$copykat_prediction, ignore.case=TRUE), na.rm=TRUE),
      predicted_diploid=sum(grepl("diploid", pred$copykat_prediction, ignore.case=TRUE), na.rm=TRUE), stringsAsFactors=FALSE)
  }
  rm(p_obj,r_obj,p,r,mat,md,ck); gc(verbose=FALSE)
}
if(length(all_calls)) fwrite(rbindlist(all_calls, fill=TRUE), file.path(out_dir,"PDAC_copykat_all_cell_calls.tsv"), sep="\t", na="NA")
if(length(run_meta)) fwrite(rbindlist(run_meta, fill=TRUE), file.path(out_dir,"PDAC_copykat_run_summary.tsv"), sep="\t", na="NA")
writeLines(capture.output(sessionInfo()), file.path(out_dir,"copykat_sessionInfo.txt"))
