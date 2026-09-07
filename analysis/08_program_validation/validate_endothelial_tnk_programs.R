options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR)})
root <- "D:/PDAC_P1"
in_cand <- file.path(root,"analysis/07_candidate_screen")
pb_dir <- file.path(root,"analysis/05_pseudobulk")
out_dir <- file.path(root,"analysis/08_program_validation")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

cand <- fread(file.path(in_cand,"candidate_genes_stable_high_confidence_all.tsv"))[candidate_status=="screening_candidate"]
programs <- split(cand$gene_id,cand$major_lineage)
programs <- programs[names(programs) %in% c("endothelial","T_NK")]
marker_sets <- list(endothelial=c("PECAM1","VWF","EMCN"),T_NK=c("CD3D","NKG7","CD8A","GNLY"))
score_mat <- function(counts,genes) {
  y <- DGEList(counts=counts); z <- cpm(y,log=TRUE,prior.count=2)
  genes <- intersect(genes,rownames(z)); if(length(genes)<3) return(rep(NA_real_,ncol(z)))
  zz <- t(scale(t(z[genes,,drop=FALSE]))); colMeans(zz,na.rm=TRUE)
}
sc_rows <- list(); sc_sum <- list()
for(ds in c("GSE155698","GSE212966")) {
  counts <- readRDS(file.path(pb_dir,paste0(ds,"_patient_pseudobulk_counts.rds")))
  md <- fread(file.path(pb_dir,paste0(ds,"_patient_pseudobulk_metadata.tsv")))
  for(lineage in names(programs)) {
    sub <- md[major_lineage==lineage & review_flag!="exclude_or_manual_review"]
    ids <- intersect(sub$group_id,colnames(counts)); sub <- sub[group_id %in% ids]
    if(!nrow(sub)) next
    s <- score_mat(counts[,sub$group_id,drop=FALSE],programs[[lineage]])
    one <- copy(sub); one[,`:=`(program=lineage,program_score=s)]; sc_rows <- append(sc_rows,list(one))
    for(ctx in unique(one$tissue_context)) {
      a <- one[tissue_context==ctx,program_score]
      sc_sum <- append(sc_sum,list(data.frame(dataset=ds,program=lineage,context=ctx,n=length(a),median_score=median(a,na.rm=TRUE),stringsAsFactors=FALSE)))
    }
    if(all(c("primary_tumor","adjacent_normal") %in% one$tissue_context)) {
      tum <- one[tissue_context=="primary_tumor",program_score]; ref <- one[tissue_context=="adjacent_normal",program_score]
      wt <- if(length(tum)>=2 && length(ref)>=2) suppressWarnings(wilcox.test(tum,ref,exact=FALSE)) else NULL
      sc_sum <- append(sc_sum,list(data.frame(dataset=ds,program=lineage,context="tumor_vs_adjacent",n_tumor=length(tum),n_adjacent=length(ref),median_tumor=median(tum,na.rm=TRUE),median_adjacent=median(ref,na.rm=TRUE),delta=median(tum,na.rm=TRUE)-median(ref,na.rm=TRUE),wilcox_P=if(is.null(wt)) NA_real_ else wt$p.value,stringsAsFactors=FALSE)))
    }
  }
}
if(length(sc_rows)) fwrite(rbindlist(sc_rows,fill=TRUE),file.path(out_dir,"scRNA_patient_program_scores.tsv"),sep="\t",na="NA")
if(length(sc_sum)) fwrite(rbindlist(sc_sum,fill=TRUE),file.path(out_dir,"scRNA_program_score_summary.tsv"),sep="\t",na="NA")

tcga_counts <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_primary_tumor_counts_unstranded.rds"))
ann <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_gene_annotation.rds"))
symbols <- ann$gene_name[match(rownames(tcga_counts),ann$gene_id)]; keep <- !is.na(symbols) & symbols!=""; tcga_counts <- tcga_counts[keep,,drop=FALSE]; rownames(tcga_counts) <- make.unique(symbols[keep])
tcga_rows <- list(); tcga_cor <- list()
for(lineage in names(programs)) {
  s <- score_mat(tcga_counts,programs[[lineage]])
  one <- data.frame(sample_id=colnames(tcga_counts),program=lineage,program_score=s,stringsAsFactors=FALSE)
  tcga_rows <- append(tcga_rows,list(one))
  for(mk in marker_sets[[lineage]]) {
    ix <- which(rownames(tcga_counts)==mk)[1]
    if(!is.na(ix) && length(ix)) {
      y <- cpm(DGEList(counts=tcga_counts),log=TRUE,prior.count=2)[ix,]
      tcga_cor <- append(tcga_cor,list(data.frame(program=lineage,marker=mk,n=length(s),spearman_rho=cor(s,y,use="pairwise.complete.obs",method="spearman"),spearman_P=cor.test(s,y,method="spearman",exact=FALSE)$p.value,stringsAsFactors=FALSE)))
    }
  }
}
if(length(tcga_rows)) fwrite(rbindlist(tcga_rows),file.path(out_dir,"TCGA_program_scores.tsv"),sep="\t",na="NA")
if(length(tcga_cor)) fwrite(rbindlist(tcga_cor),file.path(out_dir,"TCGA_program_marker_correlations.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"program_validation_sessionInfo.txt"))
