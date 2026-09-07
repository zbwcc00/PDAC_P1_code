options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
root <- "D:/PDAC_P1"
qc_dir <- file.path(root, "analysis/02_scrna_qc/seurat_checkpoints")
ann_dir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")
out_dir <- file.path(root, "analysis/03_scrna_annotation")
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
marker_sets <- list(
  epithelial=c("EPCAM","KRT8","KRT18","KRT19","KRT7","KRT17","MUC1"),
  malignant_epithelial=c("EPCAM","KRT8","KRT18","KRT19","KRT7","KRT17","MSLN","CEACAM6","MUC1","S100A10"),
  fibroblast_CAF=c("COL1A1","COL1A2","DCN","LUM","COL3A1","SPARC","FAP","PDGFRA","PDGFRB","ACTA2"),
  myeloid=c("LST1","TYROBP","FCER1G","CTSS","AIF1","LILRB1","LGALS3","LYZ","S100A8","S100A9"),
  T_NK=c("CD3D","CD3E","TRBC1","TRBC2","NKG7","GNLY","KLRD1"),
  CD8=c("CD3D","CD3E","CD8A","CD8B","NKG7","GZMK","GZMB","PRF1"),
  B_cell=c("CD79A","MS4A1","CD37","CD74","HLA-DRA","CD79B"),
  endothelial=c("EMCN","PECAM1","VWF","KDR","RAMP2","ESAM","ENG"),
  mast=c("KIT","TPSAB1","MS4A2","CPA3"),
  acinar=c("PRSS1","REG1A","REG1B","AMY2A","CPA1","CTRB1"))
module_sets <- list(malignant_epithelial_program=marker_sets$malignant_epithelial,
  CAF_program=marker_sets$fibroblast_CAF, myeloid_program=marker_sets$myeloid,
  CD8_activation=c("CD3D","CD3E","CD8A","CD8B","NKG7","GNLY","GZMK","GZMB","PRF1"),
  CD8_exhaustion=c("PDCD1","LAG3","TIGIT","TOX","HAVCR2","CTLA4","TNFRSF9","ENTPD1"))
score_genes <- function(expr, genes) { genes <- intersect(genes, rownames(expr)); if(!length(genes)) return(rep(0,ncol(expr))); Matrix::colMeans(expr[genes,,drop=FALSE]) }
classify <- function(scores) { best <- max.col(scores,ties.method="first"); bs <- scores[cbind(seq_len(nrow(scores)),best)]; s <- t(apply(scores,1,sort,decreasing=TRUE)); lab <- colnames(scores)[best]; lab[bs<0.25 | s[,1]-s[,2]<0.08] <- "unknown_ambiguous"; lab }
annotate_one <- function(infile, outfile) {
  obj <- readRDS(infile); if("RNA" %in% names(obj@assays)) obj <- JoinLayers(obj, assay="RNA"); obj <- NormalizeData(obj, verbose=FALSE)
  expr <- GetAssayData(obj, assay="RNA", layer="data"); ls <- as.matrix(sapply(marker_sets, function(g) score_genes(expr,g))); colnames(ls) <- names(marker_sets); obj$major_lineage <- classify(ls)
  for(j in seq_len(ncol(ls))) obj[[paste0("score_",colnames(ls)[j])]] <- ls[,j]
  for(nm in names(module_sets)) obj[[paste0("module_",nm)]] <- score_genes(expr,module_sets[[nm]])
  saveRDS(obj,outfile,compress=FALSE); invisible(obj)
}
for(key in c("GSE155698__primary_tumor","GSE212966__primary_tumor","GSE212966__adjacent_normal")) {
  out <- file.path(ann_dir,paste0(key,".rds")); if(!file.exists(out) || file.info(out)$size==0) { message("Annotating ",key); annotate_one(file.path(qc_dir,paste0(key,".rds")),out) }
}

files <- list.files(ann_dir, pattern="\\.rds$", full.names=TRUE)
cell_rows <- list()
for(f in files) {
  obj <- tryCatch(readRDS(f), error=function(e) NULL); if(is.null(obj)) next
  md <- as.data.frame(obj[[]]); md$cell <- rownames(md); setDT(md)
  needed <- c("cell","dataset","gsm","tissue_context","patient_id","major_lineage","module_malignant_epithelial_program","module_CAF_program","module_myeloid_program","module_CD8_activation","module_CD8_exhaustion")
  if(!all(needed %in% names(md))) next
  md$pass_singlet <- md$passes_core_QC %in% TRUE & (is.na(md$doublet_class) | md$doublet_class != "doublet")
  md <- md[pass_singlet == TRUE, ..needed]
  cell_rows <- append(cell_rows,list(md)); rm(obj); gc(verbose=FALSE)
}
cells <- rbindlist(cell_rows, fill=TRUE)
fwrite(cells,file.path(out_dir,"PDAC_scRNA_cell_major_lineage_scores.tsv"),sep="\t",na="NA")
denom <- cells[,.(total_singlet_cells=.N),by=.(dataset,gsm,tissue_context,patient_id)]
counts <- dcast(cells,dataset+gsm+tissue_context+patient_id~major_lineage,fun.aggregate=length,value.var="cell",fill=0)
wide <- merge(denom,counts,by=c("dataset","gsm","tissue_context","patient_id"),all=TRUE)
lineage_cols <- setdiff(names(wide),c("dataset","gsm","tissue_context","patient_id","total_singlet_cells"))
for(nm in lineage_cols) wide[[paste0("prop_",nm)]] <- wide[[nm]]/wide$total_singlet_cells
mods <- cells[,.(malignant_epithelial_program=mean(module_malignant_epithelial_program),CAF_program=mean(module_CAF_program),myeloid_program=mean(module_myeloid_program),CD8_activation=mean(module_CD8_activation),CD8_exhaustion=mean(module_CD8_exhaustion)),by=.(dataset,gsm,tissue_context,patient_id)]
wide <- merge(wide,mods,by=c("dataset","gsm","tissue_context","patient_id"),all=TRUE)
fwrite(wide,file.path(out_dir,"PDAC_scRNA_patient_lineage_summary.tsv"),sep="\t",na="NA")
primary <- wide[tissue_context=="primary_tumor"]; cor_rows <- list()
safe_cor <- function(a,b) if(length(unique(a))>2 && length(unique(b))>2) suppressWarnings(cor(a,b,method="spearman",use="pairwise.complete.obs")) else NA_real_
for(ds in unique(primary$dataset)) {x<-primary[dataset==ds]; if(nrow(x)>=6) cor_rows<-append(cor_rows,list(data.frame(dataset=ds,n_patients=nrow(x),rho_malignant_CAF=safe_cor(x$malignant_epithelial_program,x$CAF_program),rho_malignant_myeloid=safe_cor(x$malignant_epithelial_program,x$myeloid_program),rho_malignant_CD8_exhaustion=safe_cor(x$malignant_epithelial_program,x$CD8_exhaustion),rho_CAF_CD8_activation=safe_cor(x$CAF_program,x$CD8_activation)))) }
if(length(cor_rows)) fwrite(rbindlist(cor_rows),file.path(out_dir,"PDAC_scRNA_primary_ecology_correlations.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"annotation_sessionInfo.txt"))
