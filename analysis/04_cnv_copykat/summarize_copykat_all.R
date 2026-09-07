options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1/analysis/04_cnv_copykat"
datasets <- c("GSE154778","GSE155698","GSE212966")
all <- list(); run <- list()
for(ds in datasets) {
  dir <- file.path(root,ds); inp <- tryCatch(readRDS(file.path(dir,"copykat_input_sampled.rds")),error=function(e) NULL)
  pred_file <- file.path(dir,paste0("PDAC_",ds,"_copykat_prediction.txt"))
  if(is.null(inp) || !file.exists(pred_file)) {run<-append(run,list(data.frame(dataset=ds,status="missing_prediction")));next}
  pred <- fread(pred_file); setnames(pred,c("cell.names","copykat.pred"),c("cell","copykat_prediction")); md<-as.data.table(inp$metadata)
  calls<-merge(pred,md[,.(cell,dataset,gsm,tissue_context,patient_id,cnv_role)],by="cell",all.x=TRUE)
  if(ds=="GSE155698") calls[,patient_id:=sub("[A-Z]$","",patient_id)]
  calls[,predicted_aneuploid:=tolower(copykat_prediction)=="aneuploid"]
  fwrite(calls,file.path(dir,"copykat_cell_calls.tsv"),sep="\t",na="NA")
  sm<-calls[,.(sampled_cells=.N,aneuploid_cells=sum(predicted_aneuploid,na.rm=TRUE),aneuploid_fraction=mean(predicted_aneuploid,na.rm=TRUE)),by=.(dataset,gsm,tissue_context,patient_id,cnv_role)]
  fwrite(sm,file.path(dir,"copykat_sample_summary.tsv"),sep="\t",na="NA"); all<-append(all,list(sm)); run<-append(run,list(data.frame(dataset=ds,status="core_prediction_complete",sampled_cells=nrow(calls),aneuploid_fraction=mean(calls$predicted_aneuploid,na.rm=TRUE))))
}
if(length(all)) fwrite(rbindlist(all,fill=TRUE),file.path(root,"PDAC_copykat_all_sample_summary.tsv"),sep="\t",na="NA")
fwrite(rbindlist(run,fill=TRUE),file.path(root,"PDAC_copykat_run_summary.tsv"),sep="\t",na="NA")
