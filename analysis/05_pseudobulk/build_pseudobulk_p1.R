options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
root <- "D:/PDAC_P1"
qc_dir <- file.path(root,"analysis/02_scrna_qc/seurat_checkpoints")
ann_dir <- file.path(root,"analysis/03_scrna_annotation")
out_dir <- file.path(root,"analysis/05_pseudobulk")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
ann <- fread(file.path(ann_dir,"PDAC_scRNA_cell_major_lineage_scores.tsv"))
flags <- fread(file.path(root,"analysis/02_scrna_qc/PDAC_scRNA_sample_review_flags.tsv"))[,.(dataset,gsm,review_flag)]
datasets <- list(
 GSE154778=c("GSE154778__primary_tumor.rds","GSE154778__metastasis.rds"),
 GSE155698=c("GSE155698__primary_tumor.rds","GSE155698__adjacent_normal.rds"),
 GSE212966=c("GSE212966__primary_tumor.rds","GSE212966__adjacent_normal.rds"))
all_meta <- list()
for(ds in names(datasets)) {
  message("Building pseudobulk ",ds); groups <- list(); meta_rows <- list(); gene_universe <- NULL
  for(fn in datasets[[ds]]) {
    obj <- readRDS(file.path(qc_dir,fn)); if("RNA" %in% names(obj@assays)) obj <- JoinLayers(obj,assay="RNA")
    md <- as.data.frame(obj[[]]); md$cell <- rownames(md); setDT(md)
    md <- md[passes_core_QC %in% TRUE & (is.na(doublet_class) | doublet_class != "doublet")]
    context <- unique(md$tissue_context); context <- context[!is.na(context)][1]
    keep_ann <- ann[dataset==ds & gsm %in% unique(md$gsm), .(cell,major_lineage,gsm,tissue_context,patient_id)]
    join <- merge(md[,.(cell,gsm,tissue_context,patient_id)],keep_ann,by=c("cell","gsm","tissue_context","patient_id"))
    join <- join[major_lineage != "unknown_ambiguous"]
    mat <- LayerData(obj,assay="RNA",layer="counts")
    gene_universe <- if(is.null(gene_universe)) rownames(mat) else intersect(gene_universe, rownames(mat))
    for(g in unique(join[,.(gsm,major_lineage)]$major_lineage)) {
      for(s in unique(join[major_lineage==g]$gsm)) {
        ids <- join[major_lineage==g & gsm==s]$cell; ids <- intersect(ids,colnames(mat)); if(length(ids)<20) next
        key <- paste(ds,s,g,sep="__"); groups[[key]] <- Matrix::rowSums(mat[,ids,drop=FALSE])
        one <- join[major_lineage==g & gsm==s][1]
        rf <- flags[dataset==ds & gsm==s]$review_flag; if(!length(rf)) rf <- "acceptable"
        meta_rows <- append(meta_rows,list(data.frame(group_id=key,dataset=ds,gsm=s,tissue_context=one$tissue_context,patient_id=one$patient_id,major_lineage=g,n_cells=length(ids),review_flag=rf[1],analysis_eligible=one$tissue_context=="primary_tumor" && rf[1]!="exclude_or_manual_review",stringsAsFactors=FALSE)))
      }
    }
    rm(obj,mat,md,join); gc(verbose=FALSE)
  }
  if(!length(groups)) next
  genes <- intersect(gene_universe,names(groups[[1]])); if(length(genes)<5000) genes <- names(groups[[1]])
  mat_pb <- do.call(cbind,lapply(groups,function(v) {z<-v[genes]; as.numeric(z)})); rownames(mat_pb)<-genes; colnames(mat_pb)<-names(groups)
  meta_pb <- rbindlist(meta_rows,fill=TRUE)
  saveRDS(mat_pb,file.path(out_dir,paste0(ds,"_pseudobulk_counts.rds")),compress="xz")
  fwrite(meta_pb,file.path(out_dir,paste0(ds,"_pseudobulk_metadata.tsv")),sep="\t",na="NA")
  all_meta <- append(all_meta,list(meta_pb)); rm(mat_pb,groups,meta_pb); gc(verbose=FALSE)
}
if(length(all_meta)) fwrite(rbindlist(all_meta,fill=TRUE),file.path(out_dir,"PDAC_pseudobulk_metadata_all.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"pseudobulk_sessionInfo.txt"))
