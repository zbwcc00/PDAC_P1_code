options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(copykat)})
set.seed(20260902)
root <- "D:/PDAC_P1"
qc_dir <- file.path(root,"analysis/02_scrna_qc/seurat_checkpoints")
ann_file <- file.path(root,"analysis/03_scrna_annotation/PDAC_scRNA_cell_major_lineage_scores.tsv")
out_root <- file.path(root,"analysis/04_cnv_copykat_epithelial"); dir.create(out_root,recursive=TRUE,showWarnings=FALSE)
ann <- fread(ann_file)[major_lineage %in% c("epithelial","malignant_epithelial"),.(cell,dataset,gsm,tissue_context,patient_id,major_lineage)]
cfg <- list(GSE154778=list(primary="GSE154778__primary_tumor.rds",reference="GSE154778__metastasis.rds"),GSE155698=list(primary="GSE155698__primary_tumor.rds",reference="GSE155698__adjacent_normal.rds"),GSE212966=list(primary="GSE212966__primary_tumor.rds",reference="GSE212966__adjacent_normal.rds"))
extract <- function(obj, dataset_name, max_n, role) {
  if("RNA"%in%names(obj@assays)) obj<-JoinLayers(obj,assay="RNA")
  md<-as.data.frame(obj[[]]);md$cell<-rownames(md);setDT(md);md<-md[passes_core_QC%in%TRUE & (is.na(doublet_class)|doublet_class!="doublet")]
  j<-merge(md[,.(cell,gsm,tissue_context,patient_id)],ann[dataset==dataset_name],by=c("cell","gsm","tissue_context","patient_id"))
  if(!nrow(j)) return(list(counts=NULL,metadata=NULL))
  j[,cnv_role:=role]; j<-j[, .SD[sample.int(.N,min(.N,max_n))], by=gsm]
  mat<-LayerData(obj,assay="RNA",layer="counts")[,j$cell,drop=FALSE]; colnames(mat)<-j$cell
  list(counts=mat,metadata=j)
}
run_meta<-list()
for(ds in names(cfg)) {
  message("Epithelial CopyKAT ",ds); po<-readRDS(file.path(qc_dir,cfg[[ds]]$primary)); ro<-readRDS(file.path(qc_dir,cfg[[ds]]$reference))
  p<-extract(po,ds,200,"tumor_or_metastasis"); r<-extract(ro,ds,100,"reference_context")
  if(is.null(p$counts)||is.null(r$counts)){run_meta<-append(run_meta,list(data.frame(dataset=ds,status="no_epithelial_cells")));next}
  genes<-intersect(rownames(p$counts),rownames(r$counts)); p$counts<-p$counts[genes,,drop=FALSE];r$counts<-r$counts[genes,,drop=FALSE]
  mat<-cbind(p$counts,r$counts); keep<-Matrix::rowSums(mat>0)>=10 & !grepl("^(MT-|RPL|RPS|MTRNR)",rownames(mat),ignore.case=TRUE); mat<-mat[keep,,drop=FALSE]
  md<-rbind(p$metadata,r$metadata); rownames(md)<-md$cell; run_dir<-file.path(out_root,ds);dir.create(run_dir,recursive=TRUE,showWarnings=FALSE);saveRDS(list(counts=mat,metadata=md),file.path(run_dir,"copykat_input_sampled.rds"),compress="xz")
  old<-getwd();setwd(run_dir)
  ck<-tryCatch(copykat(rawmat=as.matrix(mat),id.type="S",cell.line="no",ngene.chr=5,min.gene.per.cell=200,LOW.DR=0.05,UP.DR=0.1,win.size=25,sam.name=paste0("PDAC_",ds,"_epithelial"),output.seg="TRUE",plot.genes="FALSE",genome="hg20",n.cores=1),error=function(e)e)
  setwd(old)
  if(inherits(ck,"error")){run_meta<-append(run_meta,list(data.frame(dataset=ds,status="error",error=conditionMessage(ck))));next}
  pred_file<-file.path(run_dir,paste0("PDAC_",ds,"_epithelial_copykat_prediction.txt")); if(file.exists(pred_file)){pred<-fread(pred_file);setnames(pred,c("cell.names","copykat.pred"),c("cell","copykat_prediction"));calls<-merge(pred,md[,.(cell,dataset,gsm,tissue_context,patient_id,cnv_role,major_lineage)],by="cell",all.x=TRUE);calls[,predicted_aneuploid:=tolower(copykat_prediction)=="aneuploid"];fwrite(calls,file.path(run_dir,"copykat_cell_calls.tsv"),sep="\t",na="NA");sm<-calls[,.(sampled_cells=.N,aneuploid_cells=sum(predicted_aneuploid,na.rm=TRUE),aneuploid_fraction=mean(predicted_aneuploid,na.rm=TRUE)),by=.(dataset,gsm,tissue_context,patient_id,cnv_role)];fwrite(sm,file.path(run_dir,"copykat_sample_summary.tsv"),sep="\t",na="NA");run_meta<-append(run_meta,list(data.frame(dataset=ds,status="core_prediction_complete",sampled_cells=nrow(calls),aneuploid_fraction=mean(calls$predicted_aneuploid,na.rm=TRUE))))}
  rm(po,ro,p,r,mat,md,ck);gc(verbose=FALSE)
}
if(length(run_meta)) fwrite(rbindlist(run_meta,fill=TRUE),file.path(out_root,"PDAC_copykat_epithelial_run_summary.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out_root,"copykat_epithelial_sessionInfo.txt"))
