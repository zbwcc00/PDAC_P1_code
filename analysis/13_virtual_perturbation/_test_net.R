suppressPackageStartupMessages({library(Seurat);library(scTenifoldNet);library(Matrix)})
o <- readRDS('D:/PDAC_P1/analysis/03_scrna_annotation/seurat_annotated/GSE154778__primary_tumor.rds')
sel <- rownames(o@meta.data)[o@meta.data$major_lineage=='endothelial']; mat <- GetAssayData(o,assay='RNA',layer='counts')[,sel,drop=FALSE]
set.seed(1); v<-apply(as.matrix(mat),1,var); keep<-names(sort(v,decreasing=TRUE))[seq_len(min(100,length(v)))]; mat<-mat[unique(c(keep,'EPAS1')), ,drop=FALSE]; Y<-mat; Y['EPAS1',]<-max(1,as.numeric(quantile(mat['EPAS1',],.95))*4); z<-scTenifoldNet(mat,Y,qc=FALSE,nc_nNet=2,nc_nCells=50,nc_nComp=2,ma_nDim=2,nCores=1); print(names(z)); print(head(z$diffRegulation)); print(str(z$diffRegulation))
