options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
root <- 'D:/PDAC_P1'
outdir <- file.path(root,'analysis/13_virtual_perturbation')
resdir <- file.path(outdir,'results')
objdir <- file.path(root,'analysis/03_scrna_annotation/seurat_annotated')
programs <- list(
  endothelial=c('CLDN5','KDR','PECAM1','VWF','EMCN','ESAM','RAMP2','PLVAP','ENG','EPAS1'),
  CAF_ECM=c('COL1A1','COL1A2','COL3A1','FN1','CTHRC1','DCN','LUM','SPARC'),
  myeloid_inflammation=c('S100A8','S100A9','LYZ','FCN1','CTSS','TYROBP','LST1','IL1B','CCL2','CCL3'),
  cytotoxic_TNK=c('NKG7','GNLY','GZMB','PRF1','CTSW','GZMH','CD8A'),
  T_exhaustion=c('PDCD1','CTLA4','LAG3','TIGIT','HAVCR2','TOX','ENTPD1'),
  TNF_NFkB=c('TNF','NFKB1','RELA','RELB','NFKBIA','TNFAIP3','CXCL8','ICAM1'),
  angiogenesis=c('VEGFA','KDR','PECAM1','VWF','EMCN','EPAS1','ENG'),
  TGFb_fibrosis=c('TGFB1','TGFBR1','TGFBR2','SMAD2','SMAD3','SMAD4','SERPINE1','CTHRC1','COL1A1','COL3A1')
)
summarize_network <- function(dr, dataset, celltype, target, perturbation) {
  do.call(rbind,lapply(names(programs),function(p){
    ix <- which(dr$gene %in% programs[[p]])
    if(!length(ix)) return(data.frame(dataset=dataset,cell_type=celltype,target=target,perturbation=perturbation,program=p,n_program_genes=0,n_sig_fdr=NA,mean_abs_Z=NA,mean_FC=NA,network_score=NA))
    pa <- dr$p.adj[ix]; z <- dr$Z[ix]; fc <- dr$FC[ix]
    data.frame(dataset=dataset,cell_type=celltype,target=target,perturbation=perturbation,program=p,n_program_genes=length(ix),n_sig_fdr=sum(pa<.05,na.rm=TRUE),mean_abs_Z=mean(abs(z),na.rm=TRUE),mean_FC=mean(fc,na.rm=TRUE),network_score=sum(-log10(pmax(pa,1e-16)),na.rm=TRUE))
  }))
}
proxy_delta <- function(mat,target,oe) {
  lib <- Matrix::colSums(mat); logx <- log1p(t(t(mat)/pmax(lib,1))*1e4)
  before <- sapply(programs,function(g){ix<-intersect(g,rownames(logx));if(length(ix)<2) return(NA_real_); z<-t(scale(t(as.matrix(logx[ix,,drop=FALSE]))));median(colMeans(z,na.rm=TRUE),na.rm=TRUE)})
  pert <- mat
  if(oe) { x<-as.numeric(mat[target,]); nz<-x[x>0]; level<-if(length(nz)) max(1,as.numeric(quantile(nz,.95))*4) else 1; pert[target,]<-pmax(x,level) } else pert[target,]<-0
  lib2 <- Matrix::colSums(pert); logp <- log1p(t(t(pert)/pmax(lib2,1))*1e4)
  after <- sapply(programs,function(g){ix<-intersect(g,rownames(logp));if(length(ix)<2) return(NA_real_); z<-t(scale(t(as.matrix(logp[ix,,drop=FALSE]))));median(colMeans(z,na.rm=TRUE),na.rm=TRUE)})
  data.frame(program=names(programs),before=as.numeric(before),after=as.numeric(after),delta=as.numeric(after-before),proxy=ifelse(oe,'OE','KO'))
}
files <- list.files(resdir,pattern='^GSE.*__(KO|OE)\\.rds$',full.names=TRUE)
net_rows <- lapply(files,function(f){bn<-basename(f); z<-strsplit(sub('\\.rds$','',bn),'__',fixed=TRUE)[[1]]; r<-readRDS(f); x<-summarize_network(r$diffRegulation,z[1],z[2],z[3],z[4]); rm(r);gc();x})
net <- do.call(rbind,net_rows)
write.table(net,file.path(outdir,'virtual_perturbation_network_program_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
ds <- c(GSE154778_primary='GSE154778__primary_tumor.rds',GSE155698_primary='GSE155698__primary_tumor.rds',GSE212966_primary='GSE212966__primary_tumor.rds')
targets <- c('EPAS1','TACC1','MARCKS','HERPUD1'); ctys <- c('endothelial','fibroblast_CAF','myeloid','T_NK')
proxy_rows <- list(); cov_rows <- list(); run_rows <- list()
for(dn in names(ds)) {
  obj <- readRDS(file.path(objdir,ds[[dn]])); assay <- if('RNA' %in% names(obj@assays)) 'RNA' else DefaultAssay(obj); counts <- GetAssayData(obj,assay=assay,layer='counts'); ct <- as.character(obj@meta.data$major_lineage); names(ct)<-rownames(obj@meta.data)
  for(cty in ctys) {
    cells <- names(ct)[ct==cty]; sub <- counts[,cells,drop=FALSE]; sub <- sub[,Matrix::colSums(sub)>=100,drop=FALSE]
    for(target in targets) {
      if(!(target %in% rownames(sub))) {run_rows[[length(run_rows)+1L]]<-data.frame(dataset=dn,cell_type=cty,target=target,perturbation='SKIP',status='gene_missing',n_cells=ncol(sub),n_genes=NA);next}
      det <- sum(as.numeric(sub[target,])>0); cov_rows[[length(cov_rows)+1L]]<-data.frame(dataset=dn,cell_type=cty,target=target,n_cells=ncol(sub),n_detected=det,detection_rate=det/max(1,ncol(sub)),n_genes=nrow(sub))
      for(pt in c('KO','OE')) {
        key <- paste(dn,cty,target,pt,sep='__'); ok <- file.exists(file.path(resdir,paste0(key,'.rds')))
        run_rows[[length(run_rows)+1L]]<-data.frame(dataset=dn,cell_type=cty,target=target,perturbation=pt,status=ifelse(ok,'OK','MISSING_RDS'),n_cells=ncol(sub),n_genes=NA)
        if(ok && det>=5) { keep <- intersect(unique(c(unlist(programs),target)),rownames(sub)); pr <- proxy_delta(sub[keep,,drop=FALSE],target,pt=='OE'); pr$dataset<-dn;pr$cell_type<-cty;pr$target<-target;pr$n_cells<-ncol(sub);pr$n_genes<-length(keep);proxy_rows[[length(proxy_rows)+1L]]<-pr }
      }
    }
  }
  rm(obj,counts);gc()
}
proxy <- do.call(rbind,proxy_rows); cov <- do.call(rbind,cov_rows); runs <- do.call(rbind,run_rows)
write.table(proxy,file.path(outdir,'virtual_perturbation_proxy_program_delta.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(cov,file.path(outdir,'virtual_perturbation_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(runs,file.path(outdir,'virtual_perturbation_run_manifest.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
agg <- aggregate(cbind(network_score,mean_abs_Z,n_sig_fdr)~target+perturbation+cell_type+program,data=net,FUN=function(x) c(n=sum(!is.na(x)),median=median(x,na.rm=TRUE),mean=mean(x,na.rm=TRUE)))
write.table(agg,file.path(outdir,'virtual_perturbation_network_aggregate.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'summary_sessionInfo.txt'))
