options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Matrix);library(data.table)})
source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- file.path(pdac_paths()$analysis_root, "05_pseudobulk")
datasets <- c("GSE154778","GSE155698","GSE212966")
all_meta <- list()
for(ds in datasets) {
  mat <- readRDS(file.path(root,paste0(ds,"_pseudobulk_counts.rds")))
  md <- fread(file.path(root,paste0(ds,"_pseudobulk_metadata.tsv")))
  md[dataset=="GSE155698",patient_id:=sub("[A-Z]$","",patient_id)]
  md <- md[group_id %in% colnames(mat)]
  groups <- list(); rows <- list()
  for(g in unique(md$major_lineage)) for(ctx in unique(md[major_lineage==g]$tissue_context)) for(pid in unique(md[major_lineage==g & tissue_context==ctx]$patient_id)) {
    x <- md[major_lineage==g & tissue_context==ctx & patient_id==pid]
    if(sum(x$n_cells)<20) next
    ids <- intersect(x$group_id,colnames(mat)); if(!length(ids)) next
    key <- paste(ds,pid,ctx,g,sep="__"); groups[[key]] <- Matrix::rowSums(mat[,ids,drop=FALSE])
    rows <- append(rows,list(data.frame(group_id=key,dataset=ds,patient_id=pid,tissue_context=ctx,major_lineage=g,gsm_list=paste(x$gsm,collapse=";"),n_cells=sum(x$n_cells),review_flag=paste(unique(x$review_flag),collapse=";"),analysis_eligible=all(x$analysis_eligible),stringsAsFactors=FALSE)))
  }
  if(!length(groups)) next
  genes <- rownames(mat); pb <- do.call(cbind,lapply(groups,function(v) as.numeric(v[genes]))); rownames(pb)<-genes; colnames(pb)<-names(groups)
  out_meta <- rbindlist(rows,fill=TRUE)
  saveRDS(pb,file.path(root,paste0(ds,"_patient_pseudobulk_counts.rds")),compress="xz")
  fwrite(out_meta,file.path(root,paste0(ds,"_patient_pseudobulk_metadata.tsv")),sep="\t",na="NA")
  all_meta <- append(all_meta,list(out_meta)); rm(mat,md,pb,groups); gc(verbose=FALSE)
}
if(length(all_meta)) fwrite(rbindlist(all_meta,fill=TRUE),file.path(root,"PDAC_patient_pseudobulk_metadata_all.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(root,"patient_pseudobulk_sessionInfo.txt"))
