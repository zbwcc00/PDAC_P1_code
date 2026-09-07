options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(clusterProfiler);library(org.Hs.eg.db)})
root <- "D:/PDAC_P1/analysis/06_pseudobulk_DE"
out_root <- "D:/PDAC_P1/analysis/07_candidate_screen"
dir.create(out_root,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out_root,"candidate_genes"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out_root,"candidate_programs"),recursive=TRUE,showWarnings=FALSE)
datasets <- c("GSE155698","GSE212966")
lineages <- intersect(
  fread(file.path(root,"celltype_DE_design_audit.tsv"))[dataset==datasets[1] & status=="run",major_lineage],
  fread(file.path(root,"celltype_DE_design_audit.tsv"))[dataset==datasets[2] & status=="run",major_lineage])
gene_out <- list(); program_out <- list(); audit_out <- list()
contaminant_genes <- c("PRSS1","PRSS2","REG1A","REG1B","AMY2A","CPA1","CTRB1","CTRB2","CELA3A","CELA3B")

for(lineage in sort(unique(lineages))) {
  dat <- list(); ok <- TRUE
  for(ds in datasets) {
    safe <- gsub("[^A-Za-z0-9]+","_",lineage)
    ep <- file.path(root,"edgeR",ds,paste0(safe,"_DE.tsv"))
    dp <- file.path(root,"DESeq2",ds,paste0(safe,"_DE.tsv"))
    if(!file.exists(ep) || !file.exists(dp)) {ok <- FALSE; break}
    e <- fread(ep)[,.(gene_id,edge_logFC=logFC,edge_FDR=FDR,edge_P=PValue)]
    d <- fread(dp)[,.(gene_id,deseq_logFC=log2FoldChange,deseq_FDR=padj,deseq_P=pvalue)]
    dat[[ds]] <- merge(e,d,by="gene_id")
  }
  if(!ok) next
  m <- merge(dat[[datasets[1]]],dat[[datasets[2]]],by="gene_id",suffixes=c("_gse155698","_gse212966"))
  effect_cols <- c("edge_logFC_gse155698","deseq_logFC_gse155698","edge_logFC_gse212966","deseq_logFC_gse212966")
  fdr_cols <- c("edge_FDR_gse155698","deseq_FDR_gse155698","edge_FDR_gse212966","deseq_FDR_gse212966")
  p_cols <- c("edge_P_gse155698","deseq_P_gse155698","edge_P_gse212966","deseq_P_gse212966")
  m[,n_tests:=rowSums(!is.na(.SD)),.SDcols=effect_cols]
  m[,n_same_direction:=apply(.SD,1,function(z){z<-z[is.finite(z)&z!=0];if(length(z)<4) return(NA_integer_);as.integer(length(unique(sign(z)))==1)}),.SDcols=effect_cols]
  m[,n_FDR05:=rowSums(as.data.frame(.SD)<0.05,na.rm=TRUE),.SDcols=fdr_cols]
  m[,n_nominal05:=rowSums(as.data.frame(.SD)<0.05,na.rm=TRUE),.SDcols=p_cols]
  m[,median_logFC:=apply(.SD,1,median,na.rm=TRUE),.SDcols=effect_cols]
  m[,mean_abs_logFC:=apply(.SD,1,function(z)mean(abs(z),na.rm=TRUE)),.SDcols=effect_cols]
  m[,direction:=ifelse(median_logFC>0,"tumor_up","tumor_down")]
  m[,stable_high_confidence:=n_tests==4 & n_same_direction==1 & n_FDR05>=3 & abs(median_logFC)>=0.5]
  setorder(m,-stable_high_confidence,-n_FDR05,-mean_abs_logFC)
  m[,`:=`(major_lineage=lineage,stability_rule="same direction in 4 tests; FDR<0.05 in >=3; median |logFC|>=0.5")]
  fwrite(m,file.path(out_root,"candidate_genes",paste0(safe,"_all_merged.tsv")),sep="\t",na="NA")
  stable <- m[stable_high_confidence==TRUE]
  stable[,contamination_flag:=major_lineage!="acinar" & gene_id %in% contaminant_genes]
  stable[,candidate_status:=ifelse(contamination_flag,"exclude_contamination","screening_candidate")]
  fwrite(stable,file.path(out_root,"candidate_genes",paste0(safe,"_stable_high_confidence.tsv")),sep="\t",na="NA")
  stable_filtered <- stable[contamination_flag==FALSE]
  fwrite(stable_filtered,file.path(out_root,"candidate_genes",paste0(safe,"_stable_filtered.tsv")),sep="\t",na="NA")
  gene_out <- append(gene_out,list(stable_filtered))
  audit_out <- append(audit_out,list(data.frame(major_lineage=lineage,n_genes_merged=nrow(m),n_stable_high_confidence=nrow(stable),n_filtered=nrow(stable_filtered),n_contamination_excluded=sum(stable$contamination_flag),n_tumor_up=sum(stable_filtered$direction=="tumor_up"),n_tumor_down=sum(stable_filtered$direction=="tumor_down"),stringsAsFactors=FALSE)))
  for(dir in c("tumor_up","tumor_down")) {
    genes <- stable_filtered[direction==dir,gene_id]
    if(length(genes)<5) next
    universe <- m[is.finite(median_logFC),gene_id]
    ego <- tryCatch(enrichGO(gene=genes,universe=universe,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont="BP",pAdjustMethod="BH",qvalueCutoff=0.2,readable=TRUE),error=function(e)NULL)
    if(is.null(ego)) next
    tab <- as.data.table(as.data.frame(ego)); if(!nrow(tab)) next
    tab[,`:=`(major_lineage=lineage,direction=dir,n_stable_genes=length(genes))]
    fwrite(tab,file.path(out_root,"candidate_programs",paste0(safe,"_",dir,"_GO.tsv")),sep="\t",na="NA")
    program_out <- append(program_out,list(tab))
  }
}
if(length(gene_out)) fwrite(rbindlist(gene_out,fill=TRUE),file.path(out_root,"candidate_genes_stable_high_confidence_all.tsv"),sep="\t",na="NA")
if(length(program_out)) fwrite(rbindlist(program_out,fill=TRUE),file.path(out_root,"candidate_programs_GO_all.tsv"),sep="\t",na="NA")
if(length(audit_out)) fwrite(rbindlist(audit_out,fill=TRUE),file.path(out_root,"candidate_screen_summary.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out_root,"candidate_screen_sessionInfo.txt"))
