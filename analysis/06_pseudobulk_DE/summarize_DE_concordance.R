options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1/analysis/06_pseudobulk_DE"
audit <- fread(file.path(root,"celltype_DE_design_audit.tsv"))
run <- audit[status=="run"]
out <- list()
for(i in seq_len(nrow(run))) {
  z <- run[i]; safe <- gsub("[^A-Za-z0-9]+","_",z$major_lineage)
  ep <- file.path(root,"edgeR",z$dataset,paste0(safe,"_DE.tsv"))
  dp <- file.path(root,"DESeq2",z$dataset,paste0(safe,"_DE.tsv"))
  if(!file.exists(ep) || !file.exists(dp)) next
  e <- fread(ep); d <- fread(dp)
  e_set <- e[!is.na(FDR) & FDR<0.05, gene_id]
  d_set <- d[!is.na(padj) & padj<0.05, gene_id]
  inter <- intersect(e_set,d_set); uni <- union(e_set,d_set)
  out <- append(out,list(data.frame(dataset=z$dataset,major_lineage=z$major_lineage,comparison=z$comparison,n_reference=z$n_reference,n_tumor=z$n_tumor,edgeR_FDR05=length(e_set),DESeq2_FDR05=length(d_set),overlap_FDR05=length(inter),jaccard_FDR05=ifelse(length(uni),length(inter)/length(uni),NA_real_),stringsAsFactors=FALSE)))
}
if(length(out)) fwrite(rbindlist(out),file.path(root,"celltype_DE_concordance.tsv"),sep="\t",na="NA")
