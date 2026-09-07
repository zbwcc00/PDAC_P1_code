options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR);library(WGCNA)})
allowWGCNAThreads(nThreads=2)
root <- "D:/PDAC_P1"
out_dir <- file.path(root,"analysis/12_functional_immune/WGCNA")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
counts <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_primary_tumor_counts_unstranded.rds"))
ann <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_gene_annotation.rds"))
symbols <- ann$gene_name[match(rownames(counts),ann$gene_id)]
keep <- !is.na(symbols) & symbols!=""
counts <- counts[keep,,drop=FALSE]; rownames(counts) <- make.unique(symbols[keep])
log_expr <- cpm(DGEList(counts=counts),log=TRUE,prior.count=2)
v <- apply(log_expr,1,var); sel <- names(sort(v,decreasing=TRUE))[seq_len(min(5000,length(v)))]
dat <- t(log_expr[sel,,drop=FALSE]); g <- goodSamplesGenes(dat,verbose=0); dat <- dat[g$goodSamples, g$goodGenes, drop=FALSE]
net <- blockwiseModules(dat,power=6,TOMType="unsigned",minModuleSize=30,reassignThreshold=0,mergeCutHeight=0.25,numericLabels=TRUE,maxBlockSize=5000,randomSeed=123,verbose=0)
colors <- labels2colors(net$colors); MEs <- orderMEs(net$MEs); rownames(MEs) <- rownames(dat)
prog <- fread(file.path(root,"analysis/08_program_validation/TCGA_program_scores.tsv")); prog <- dcast(prog,sample_id~program,value.var="program_score")
gsva <- fread(file.path(root,"analysis/12_functional_immune/TCGA_immune_GSVA_scores.tsv")); gsva <- melt(gsva,id.vars="pathway",variable.name="sample_id",value.name="score"); gsva <- dcast(gsva,sample_id~pathway,value.var="score")
trait <- merge(prog,gsva,by="sample_id",all=FALSE); trait_ids <- as.character(trait[["sample_id"]]); trait$sample_id <- NULL; trait <- as.data.frame(trait); rownames(trait) <- trait_ids; cat("WGCNA names:",head(rownames(MEs),2),head(rownames(trait),2),"\n"); common <- intersect(rownames(MEs),rownames(trait)); cat("WGCNA alignment:",nrow(MEs),ncol(MEs),nrow(trait),ncol(trait),length(common),"\n"); MEs <- MEs[common,,drop=FALSE]; trait <- trait[common,,drop=FALSE]
trait_num <- as.data.frame(lapply(trait,as.numeric)); rownames(trait_num) <- rownames(trait); trait_num <- trait_num[rownames(MEs),,drop=FALSE]; trait_num <- trait_num[,vapply(trait_num,function(x)sum(is.finite(x))>3 && sd(x,na.rm=TRUE)>0,logical(1)),drop=FALSE]
if(nrow(MEs)<4 || ncol(trait_num)<1) stop("WGCNA trait alignment produced insufficient samples or traits")
cor_mat <- cor(MEs,trait_num,use="pairwise.complete.obs",method="spearman"); p_mat <- matrix(NA_real_,nrow=nrow(cor_mat),ncol=ncol(cor_mat),dimnames=dimnames(cor_mat))
for(i in seq_len(nrow(cor_mat))) for(j in seq_len(ncol(cor_mat))) {x<-MEs[,i];y<-trait_num[,j];ok<-is.finite(x)&is.finite(y);if(sum(ok)>3)p_mat[i,j]<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE)$p.value)}
cor_dt <- as.data.table(as.table(cor_mat)); setnames(cor_dt,c("module","trait","rho")); cor_dt[,p:=as.vector(p_mat)]; cor_dt[,padj:=p.adjust(p,"BH")]; fwrite(cor_dt,file.path(out_dir,"TCGA_WGCNA_module_trait_correlations.tsv"),sep="\t")
gene_dt <- data.table(gene=colnames(dat),module=colors); fwrite(gene_dt,file.path(out_dir,"TCGA_WGCNA_gene_modules.tsv"),sep="\t")
eig_dt <- as.data.table(MEs,keep.rownames="sample_id"); fwrite(eig_dt,file.path(out_dir,"TCGA_WGCNA_module_eigengenes.tsv"),sep="\t")
saveRDS(net,file.path(out_dir,"TCGA_WGCNA_network.rds")); writeLines(capture.output(sessionInfo()),file.path(out_dir,"WGCNA_sessionInfo.txt"))
