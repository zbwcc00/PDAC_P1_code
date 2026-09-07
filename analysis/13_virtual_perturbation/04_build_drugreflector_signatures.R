options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
root <- 'D:/PDAC_P1'
outdir <- file.path(root,'analysis/13_virtual_perturbation/drugreflector')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
datasets <- c(GSE154778_primary='GSE154778__primary_tumor.rds',GSE155698_primary='GSE155698__primary_tumor.rds',GSE212966_primary='GSE212966__primary_tumor.rds')
targets <- c('EPAS1','MARCKS','HERPUD1')
strata <- list(EPAS1='endothelial',MARCKS='endothelial',HERPUD1='endothelial')
program_strata <- list(EPAS1='endothelial',MARCKS='fibroblast_CAF',HERPUD1='myeloid')
calc_v <- function(dat,hi,lo) {
  x1 <- dat[,lo,drop=FALSE]; x2 <- dat[,hi,drop=FALSE]
  m1 <- Matrix::rowMeans(x1); m2 <- Matrix::rowMeans(x2)
  v1 <- Matrix::rowMeans(x1^2)-m1^2; v2 <- Matrix::rowMeans(x2^2)-m2^2
  (m2-m1)/(sqrt(pmax(v1,0)+pmax(v2,0)) + ((v1+v2)==0))
}
all_rows <- list(); meta_rows <- list(); k <- 0L
for(dn in names(datasets)) {
  obj <- readRDS(file.path(root,'analysis/03_scrna_annotation/seurat_annotated',datasets[[dn]]))
  assay <- if('RNA' %in% names(obj@assays)) 'RNA' else DefaultAssay(obj)
  dat <- GetAssayData(obj,assay=assay,layer='data')
  md <- obj@meta.data; ct <- as.character(md$major_lineage); names(ct)<-rownames(md)
  for(target in targets) {
    for(cty in unique(c(strata[[target]],program_strata[[target]]))) {
      cells <- names(ct)[ct==cty]; cells <- intersect(cells,colnames(dat));
      if(length(cells)<50 || !(target %in% rownames(dat))) next
      x <- as.numeric(dat[target,cells]); ord <- order(x); n <- length(ord); m <- max(10,floor(n*.3)); lo <- cells[ord[seq_len(m)]]; hi <- cells[ord[(n-m+1):n]]
      if(length(unique(x))<3) next
      v <- calc_v(dat[,c(lo,hi),drop=FALSE],hi,lo); k <- k+1L; sig <- paste(target,cty,dn,sep='__')
      all_rows[[k]] <- data.frame(signature=sig,gene=rownames(dat),vscore=as.numeric(v),target=target,cell_type=cty,dataset=dn)
      meta_rows[[k]] <- data.frame(signature=sig,target=target,cell_type=cty,dataset=dn,n_low=length(lo),n_high=length(hi),definition='top30pct_vs_bottom30pct_target_expression')
    }
  }
  rm(obj,dat); gc()
}
long <- do.call(rbind,all_rows); meta <- do.call(rbind,meta_rows)
write.table(long,file.path(outdir,'drugreflector_stratum_vscores.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(meta,file.path(outdir,'drugreflector_stratum_metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
comp_rows <- list(); j <- 0L
for(target in targets) {
  z <- long[long$target==target,,drop=FALSE]
  for(cty in unique(z$cell_type)) {
    q <- z[z$cell_type==cty,,drop=FALSE]
    a <- aggregate(vscore~gene,q,FUN=function(x) median(x,na.rm=TRUE)); a$signature<-paste(target,cty,'META',sep='__'); a$target<-target; a$cell_type<-cty; a$dataset<-'META'; j<-j+1L; comp_rows[[j]]<-a
  }
  q <- z[z$cell_type %in% unique(z$cell_type),,drop=FALSE]
  a <- aggregate(vscore~gene,q,FUN=function(x) median(x,na.rm=TRUE)); a$signature<-paste(target,'COMPOSITE','META',sep='__'); a$target<-target; a$cell_type<-'COMPOSITE'; a$dataset<-'META'; j<-j+1L; comp_rows[[j]]<-a
}
comp <- do.call(rbind,comp_rows); write.table(comp,file.path(outdir,'drugreflector_vscores.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'signature_sessionInfo.txt'))
