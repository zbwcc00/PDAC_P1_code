suppressPackageStartupMessages({library(Seurat);library(scTenifoldKnk);library(Matrix)})
o <- readRDS('D:/PDAC_P1/analysis/03_scrna_annotation/seurat_annotated/GSE154778__primary_tumor.rds')
sel <- rownames(o@meta.data)[o@meta.data$major_lineage=='endothelial']
mat <- GetAssayData(o,assay='RNA',layer='counts')[,sel,drop=FALSE]
set.seed(1); g <- c('EPAS1','TACC1','MARCKS','HERPUD1','CTHRC1','CLDN5','KDR','PECAM1','VWF','EMCN','ESAM','RAMP2','PLVAP','ENG')
keep <- intersect(g,rownames(mat)); mat <- mat[keep,,drop=FALSE]
print(dim(mat)); z <- scTenifoldKnk(mat,qc=FALSE,gKO='EPAS1',nc_nNet=2,nc_nCells=50,nc_nComp=2,nCores=1); print(names(z)); print(lapply(z,class)); print(str(z$diffRegulation,max.level=1))
