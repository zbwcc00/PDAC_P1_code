options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(scTenifoldKnk); library(scTenifoldNet)})

source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- pdac_paths()$project_root
outdir <- file.path(root,'analysis/13_virtual_perturbation')
objdir <- file.path(root,'analysis/03_scrna_annotation/seurat_annotated')
dir.create(file.path(outdir,'results'),recursive=TRUE,showWarnings=FALSE)
datasets <- c(GSE154778_primary='GSE154778__primary_tumor.rds', GSE155698_primary='GSE155698__primary_tumor.rds', GSE212966_primary='GSE212966__primary_tumor.rds')
targets <- c('EPAS1','TACC1','MARCKS','HERPUD1')
celltypes <- c('endothelial','fibroblast_CAF','myeloid','T_NK')
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

score_program <- function(log_expr, genes) {
  ix <- intersect(genes,rownames(log_expr))
  if(length(ix)<2) return(rep(NA_real_,ncol(log_expr)))
  z <- t(scale(t(log_expr[ix,,drop=FALSE])))
  colMeans(z,na.rm=TRUE)
}

summarize_network <- function(dr, dataset, celltype, target, perturbation, n_cells, n_genes) {
  if(is.null(dr) || !nrow(dr)) return(data.frame())
  out <- lapply(names(programs), function(p) {
    ix <- which(dr$gene %in% programs[[p]])
    if(!length(ix)) return(data.frame(dataset=dataset,cell_type=celltype,target=target,perturbation=perturbation,n_cells=n_cells,n_genes=n_genes,program=p,n_program_genes=0,n_sig_fdr=NA,mean_abs_Z=NA,mean_FC=NA,network_score=NA))
    pa <- dr$p.adj[ix]
    z <- dr$Z[ix]
    fc <- dr$FC[ix]
    data.frame(dataset=dataset,cell_type=celltype,target=target,perturbation=perturbation,n_cells=n_cells,n_genes=n_genes,program=p,n_program_genes=length(ix),n_sig_fdr=sum(pa<0.05,na.rm=TRUE),mean_abs_Z=mean(abs(z),na.rm=TRUE),mean_FC=mean(fc,na.rm=TRUE),network_score=sum(-log10(pmax(pa,1e-300)),na.rm=TRUE))
  })
  do.call(rbind,out)
}

proxy_delta <- function(mat, target, oe=FALSE) {
  lib <- Matrix::colSums(mat)
  logx <- log1p(t(t(mat)/pmax(lib,1))*1e4)
  before <- sapply(programs, function(g) mean(score_program(logx,g),na.rm=TRUE))
  pert <- mat
  if(oe) {
    x <- as.numeric(mat[target,]); nz <- x[x>0]
    level <- if(length(nz)) max(1,as.numeric(quantile(nz,.95))*4) else 1
    pert[target,] <- pmax(x,level)
  } else pert[target,] <- 0
  lib2 <- Matrix::colSums(pert)
  logp <- log1p(t(t(pert)/pmax(lib2,1))*1e4)
  after <- sapply(programs, function(g) mean(score_program(logp,g),na.rm=TRUE))
  data.frame(program=names(programs),before=as.numeric(before),after=as.numeric(after),delta=as.numeric(after-before),proxy=ifelse(oe,'OE','KO'))
}

coverage_rows <- list(); net_rows <- list(); proxy_rows <- list(); run_rows <- list(); rr <- 0L
for(dn in names(datasets)) {
  message('Reading ',dn)
  obj <- readRDS(file.path(objdir,datasets[[dn]]))
  assay <- if('RNA' %in% names(obj@assays)) 'RNA' else DefaultAssay(obj)
  counts <- GetAssayData(obj,assay=assay,layer='counts')
  md <- obj@meta.data
  ct <- as.character(md$major_lineage); names(ct) <- rownames(md)
  for(cty in celltypes) {
    cells <- names(ct)[ct==cty]
    if(length(cells)<50) { run_rows[[length(run_rows)+1L]] <- data.frame(dataset=dn,cell_type=cty,target=NA,perturbation='SKIP',status='n_cells_lt_50',n_cells=length(cells)); next }
    sub <- counts[,cells,drop=FALSE]
    lib <- Matrix::colSums(sub); sub <- sub[,lib>=100,drop=FALSE]
    if(ncol(sub)<50) { run_rows[[length(run_rows)+1L]] <- data.frame(dataset=dn,cell_type=cty,target=NA,perturbation='SKIP',status='n_cells_after_lib_filter_lt_50',n_cells=ncol(sub)); next }
    log_sub <- log1p(t(t(sub)/pmax(Matrix::colSums(sub),1))*1e4)
    vv <- apply(as.matrix(log_sub),1,var)
    base_genes <- names(sort(vv,decreasing=TRUE))[seq_len(min(120,length(vv)))]
    wanted <- unique(c(base_genes,unlist(programs),targets))
    wanted <- intersect(wanted,rownames(sub)); sub <- sub[wanted,,drop=FALSE]
    for(target in targets) {
      if(!(target %in% rownames(sub))) { run_rows[[length(run_rows)+1L]] <- data.frame(dataset=dn,cell_type=cty,target=target,perturbation='SKIP',status='gene_missing',n_cells=ncol(sub)); next }
      detected <- sum(as.numeric(sub[target,])>0)
      coverage_rows[[length(coverage_rows)+1L]] <- data.frame(dataset=dn,cell_type=cty,target=target,n_cells=ncol(sub),n_detected=detected,detection_rate=detected/ncol(sub),n_genes=nrow(sub))
      if(detected<5) { run_rows[[length(run_rows)+1L]] <- data.frame(dataset=dn,cell_type=cty,target=target,perturbation='SKIP',status='target_detected_lt_5',n_cells=ncol(sub)); next }
      ncell_net <- min(100,ncol(sub))
      for(pt in c('KO','OE')) {
        key <- paste(dn,cty,target,pt,sep='__')
        message(key)
        res <- tryCatch({
          if(pt=='KO') {
            scTenifoldKnk(sub,qc=FALSE,gKO=target,nc_nNet=2,nc_nCells=ncell_net,nc_nComp=3,td_K=3,ma_nDim=2,nCores=1)
          } else {
            oe <- sub; x <- as.numeric(oe[target,]); nz <- x[x>0]; level <- if(length(nz)) max(1,as.numeric(quantile(nz,.95))*4) else 1; oe[target,] <- pmax(x,level)
            scTenifoldNet(sub,oe,qc=FALSE,nc_nNet=2,nc_nCells=ncell_net,nc_nComp=3,td_K=3,ma_nDim=2,nCores=1)
          }
        }, error=function(e)e)
        if(inherits(res,'error')) {
          run_rows[[length(run_rows)+1L]] <- data.frame(dataset=dn,cell_type=cty,target=target,perturbation=pt,status=paste0('ERROR:',conditionMessage(res)),n_cells=ncol(sub)); next
        }
        saveRDS(res,file.path(outdir,'results',paste0(key,'.rds')),compress='xz')
        dr <- res$diffRegulation
        net_rows[[length(net_rows)+1L]] <- summarize_network(dr,dn,cty,target,pt,ncol(sub),nrow(sub))
        proxy_rows[[length(proxy_rows)+1L]] <- transform(proxy_delta(sub,target,oe=(pt=='OE')),dataset=dn,cell_type=cty,target=target,n_cells=ncol(sub),n_genes=nrow(sub))
        run_rows[[length(run_rows)+1L]] <- data.frame(dataset=dn,cell_type=cty,target=target,perturbation=pt,status='OK',n_cells=ncol(sub),n_genes=nrow(sub))
        rm(res); gc()
      }
    }
  }
  rm(obj,counts); gc()
}
cov <- if(length(coverage_rows)) do.call(rbind,coverage_rows) else data.frame()
net <- if(length(net_rows)) do.call(rbind,net_rows) else data.frame()
proxy <- if(length(proxy_rows)) do.call(rbind,proxy_rows) else data.frame()
runs <- if(length(run_rows)) do.call(rbind,run_rows) else data.frame()
write.table(cov,file.path(outdir,'virtual_perturbation_coverage.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(net,file.path(outdir,'virtual_perturbation_network_program_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(proxy,file.path(outdir,'virtual_perturbation_proxy_program_delta.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(runs,file.path(outdir,'virtual_perturbation_run_manifest.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'virtual_perturbation_sessionInfo.txt'))
