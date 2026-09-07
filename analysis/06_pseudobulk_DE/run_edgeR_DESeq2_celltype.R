options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Matrix);library(data.table);library(edgeR);library(DESeq2)})
source(file.path(pdac_script_repo_root(), "config", "paths.R"))
root <- pdac_paths()$project_root
in_dir <- file.path(root,"analysis/05_pseudobulk")
out_dir <- file.path(root,"analysis/06_pseudobulk_DE")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
datasets <- c("GSE154778","GSE155698","GSE212966")
all_summary <- list(); all_audit <- list()

for(ds in datasets) {
  counts <- readRDS(file.path(in_dir,paste0(ds,"_patient_pseudobulk_counts.rds")))
  md <- fread(file.path(in_dir,paste0(ds,"_patient_pseudobulk_metadata.tsv")))
  stopifnot(identical(colnames(counts),md$group_id))
  md <- md[review_flag != "exclude_or_manual_review"]
  keep_cols <- intersect(md$group_id,colnames(counts)); md <- md[group_id %in% keep_cols]
  counts <- counts[,md$group_id,drop=FALSE]
  if(ds=="GSE154778") {
    md[,condition:=ifelse(tissue_context=="primary_tumor","tumor","metastasis")]
    comparison <- "primary_vs_metastasis_exploratory"
    reference_label <- "metastasis"
  } else {
    md[,condition:=ifelse(tissue_context=="primary_tumor","tumor","adjacent_normal")]
    comparison <- "primary_vs_adjacent_normal"
    reference_label <- "adjacent_normal"
  }
  md[,condition:=factor(condition,levels=c(reference_label,"tumor"))]
  for(lineage in sort(unique(md$major_lineage))) {
    sub <- md[major_lineage==lineage & !is.na(condition)]
    tab <- sub[, .N, by=condition]
    if(!all(c(reference_label,"tumor") %in% tab$condition)) next
    n_ref <- sub[condition==reference_label,.N]; n_tum <- sub[condition=="tumor",.N]
    audit <- data.frame(dataset=ds,major_lineage=lineage,comparison=comparison,n_reference=n_ref,n_tumor=n_tum,review_flags=paste(sort(unique(sub$review_flag)),collapse=";"),status=ifelse(n_ref>=2 && n_tum>=2,"run","insufficient_replicates"))
    all_audit <- append(all_audit,list(audit))
    if(n_ref<2 || n_tum<2) next
    ids <- sub$group_id; y <- DGEList(counts=counts[,ids,drop=FALSE])
    y <- y[filterByExpr(y,group=sub$condition),,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    design <- model.matrix(~0+sub$condition); colnames(design) <- levels(sub$condition)
    fit <- glmQLFit(y,design,robust=TRUE)
    contrast <- rep(0,ncol(design)); names(contrast) <- colnames(design); contrast["tumor"] <- 1; contrast[reference_label] <- -1
    qlf <- glmQLFTest(fit,contrast=contrast)
    er <- as.data.table(topTags(qlf,n=Inf,sort.by="PValue")$table,keep.rownames="gene_id")
    er[,`:=`(dataset=ds,major_lineage=lineage,comparison=comparison,n_reference=n_ref,n_tumor=n_tum,method="edgeR_QLF")]
    out_base <- file.path(out_dir,"edgeR",ds); dir.create(out_base,recursive=TRUE,showWarnings=FALSE)
    fwrite(er,file.path(out_base,paste0("",gsub("[^A-Za-z0-9]+","_",lineage),"_DE.tsv")),sep="\t",na="NA")
    dds <- tryCatch({
      cts <- round(as.matrix(counts[,ids,drop=FALSE])); rownames(cts) <- rownames(counts)
      cts <- cts[rowSums(cts)>=10,,drop=FALSE]
      coldata <- data.frame(condition=sub$condition,row.names=ids)
      z <- DESeqDataSetFromMatrix(countData=cts,colData=coldata,design=~condition)
      z <- DESeq(z,quiet=TRUE); as.data.table(as.data.frame(results(z,contrast=c("condition","tumor",reference_label))),keep.rownames="gene_id")
    },error=function(e) NULL)
    if(!is.null(dds)) {
      dds[,`:=`(dataset=ds,major_lineage=lineage,comparison=comparison,n_reference=n_ref,n_tumor=n_tum,method="DESeq2")]
      out_d <- file.path(out_dir,"DESeq2",ds); dir.create(out_d,recursive=TRUE,showWarnings=FALSE)
      fwrite(dds,file.path(out_d,paste0(gsub("[^A-Za-z0-9]+","_",lineage),"_DE.tsv")),sep="\t",na="NA")
      audit$DESeq2_status <- "complete"
    } else audit$DESeq2_status <- "error"
    all_summary <- append(all_summary,list(data.frame(dataset=ds,major_lineage=lineage,comparison=comparison,n_reference=n_ref,n_tumor=n_tum,n_genes_tested=nrow(er),edgeR_FDR05=sum(er$FDR<0.05,na.rm=TRUE),DESeq2_status=audit$DESeq2_status,stringsAsFactors=FALSE)))
  }
  rm(counts,md); gc(verbose=FALSE)
}
if(length(all_audit)) fwrite(rbindlist(all_audit,fill=TRUE),file.path(out_dir,"celltype_DE_design_audit.tsv"),sep="\t",na="NA")
if(length(all_summary)) fwrite(rbindlist(all_summary,fill=TRUE),file.path(out_dir,"celltype_DE_summary.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"DE_sessionInfo.txt"))
