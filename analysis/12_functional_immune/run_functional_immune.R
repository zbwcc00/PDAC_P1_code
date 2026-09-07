options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR);library(GSVA);library(fgsea)})
root <- "D:/PDAC_P1"
out_dir <- file.path(root,"analysis/12_functional_immune")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

cand <- fread(file.path(root,"analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
cand <- cand[candidate_status=="screening_candidate"]
programs <- split(cand$gene_id,cand$major_lineage)
programs <- programs[names(programs) %in% c("endothelial","T_NK")]
programs <- lapply(programs,unique)
immune_sets <- list(
  IFN_gamma=c("IFNG","STAT1","IRF1","CXCL9","CXCL10","GZMB","IDO1"),
  TNF_NFkB=c("TNF","NFKB1","RELA","RELB","NFKBIA","TNFAIP3","CXCL8","ICAM1"),
  cytotoxic_TNK=c("NKG7","GNLY","GZMB","PRF1","CTSW","GZMH","CD8A"),
  T_exhaustion=c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","TOX","ENTPD1"),
  CAF_ECM=c("COL1A1","COL1A2","COL3A1","FN1","CTHRC1","DCN","LUM","SPARC"),
  angiogenesis=c("VEGFA","KDR","PECAM1","VWF","EMCN","EPAS1","ENG")
)
pathway_sets <- c(immune_sets,list(
  hypoxia=c("HIF1A","EPAS1","VEGFA","CA9","LDHA","SLC2A1","ADM","ANGPTL4","NDRG1"),
  endothelial_junction=c("CLDN5","KDR","PECAM1","VWF","EMCN","ESAM","RAMP2","PLVAP","ENG"),
  TGFb_fibrosis=c("TGFB1","TGFBR1","TGFBR2","SMAD2","SMAD3","SMAD4","SERPINE1","CTHRC1","COL1A1","COL3A1"),
  glycolysis=c("HK2","PFKM","ALDOA","GAPDH","ENO1","PKM","LDHA","SLC2A1","PGK1"),
  oxidative_phosphorylation=c("NDUFS1","NDUFV1","SDHA","UQCRC1","COX4I1","ATP5F1A","ATP5MC1","CYC1"),
  antigen_presentation=c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","HLA-DRA","HLA-DRB1","CD74"),
  myeloid_inflammation=c("S100A8","S100A9","LYZ","FCN1","CTSS","TYROBP","LST1","IL1B","CCL2","CCL3")
))
score_program <- function(log_expr, genes) {
  ix <- intersect(genes,rownames(log_expr))
  if(length(ix)<2) return(rep(NA_real_,ncol(log_expr)))
  z <- t(scale(t(log_expr[ix,,drop=FALSE])))
  colMeans(z,na.rm=TRUE)
}
run_gsva <- function(log_expr, sets) {
  sets <- lapply(sets,intersect,x=rownames(log_expr)); sets <- sets[lengths(sets)>=2]
  if(!length(sets)) return(NULL)
  param <- ssgseaParam(as.matrix(log_expr),sets,minSize=2,maxSize=500,normalize=TRUE)
  as.matrix(gsva(param,verbose=FALSE))
}
run_fgsea <- function(stats, label) {
  stats <- stats[is.finite(stats)]; stats <- sort(stats,decreasing=TRUE)
  stats <- stats[!duplicated(names(stats))]
  pathways <- lapply(pathway_sets,intersect,x=names(stats)); pathways <- pathways[lengths(pathways)>=3]
  z <- fgsea(pathways=pathways,stats=stats,minSize=3,maxSize=500,nperm=10000)
  z <- as.data.table(z); z[,program:=label]; z[order(padj,-abs(NES))]
}

counts <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_primary_tumor_counts_unstranded.rds"))
ann <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_gene_annotation.rds"))
symbols <- ann$gene_name[match(rownames(counts),ann$gene_id)]
keep <- !is.na(symbols) & symbols!=""
counts <- counts[keep,,drop=FALSE]; rownames(counts) <- make.unique(symbols[keep])
log_expr <- cpm(DGEList(counts=counts),log=TRUE,prior.count=2)
program_scores <- sapply(programs,score_program,log_expr=log_expr)
immune_scores <- run_gsva(log_expr,immune_sets)
gsva_scores <- run_gsva(log_expr,pathway_sets)
fwrite(data.table(sample_id=colnames(log_expr),program_scores),file.path(out_dir,"TCGA_program_scores_functional.tsv"),sep="\t")
if(!is.null(gsva_scores)) fwrite(as.data.table(gsva_scores,keep.rownames="pathway"),file.path(out_dir,"TCGA_immune_GSVA_scores.tsv"),sep="\t")

cor_rows <- list();k <- 1; gsea_rows <- list();j <- 1
for(label in names(programs)) {
  s <- program_scores[,label]
  for(pathway in rownames(immune_scores)) {
    ct <- suppressWarnings(cor.test(s,immune_scores[pathway,],method="spearman",exact=FALSE))
    cor_rows[[k]] <- data.frame(program=label,pathway=pathway,rho=unname(ct$estimate),p=ct$p.value);k <- k+1
  }
  test_p <- sapply(seq_len(nrow(log_expr)),function(i)suppressWarnings(cor.test(log_expr[i,],s,method="spearman",exact=FALSE)$p.value))
  rho <- apply(log_expr,1,function(v)suppressWarnings(cor(v,s,method="spearman",use="pairwise.complete.obs")))
  stat <- sign(rho)*qnorm(pmax(pmin(1-test_p/2,1-1e-16),1e-16))
  names(stat) <- rownames(log_expr)
  gsea_rows[[j]] <- run_fgsea(stat,label);j <- j+1
}
fwrite(rbindlist(cor_rows),file.path(out_dir,"TCGA_program_immune_GSVA_correlations.tsv"),sep="\t")
fwrite(rbindlist(gsea_rows,fill=TRUE),file.path(out_dir,"TCGA_program_correlated_Hallmark_GSEA.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"functional_immune_sessionInfo.txt"))
