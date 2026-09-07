options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(edgeR)})
root <- "D:/PDAC_P1"
out_dir <- file.path(root,"analysis/08_program_validation")
cand <- fread(file.path(root,"analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))[candidate_status=="screening_candidate"]
programs <- split(cand$gene_id,cand$major_lineage)
programs <- programs[names(programs) %in% c("endothelial","T_NK")]
files <- c(A="D:/PDAC_P1/data/03_spatial/GSE111672/raw/GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz",B="D:/PDAC_P1/data/03_spatial/GSE111672/raw/GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz")
rows <- list(); avail <- list()
for(sample_name in names(files)) {
  x <- fread(files[[sample_name]])
  genes <- x[[1]]; mat <- as.matrix(x[,-1,with=FALSE]); storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(genes); labels <- names(x)[-1]
  lib <- colSums(mat); z <- log1p(t(t(mat)/pmax(lib,1)*1e6))
  for(program in names(programs)) {
    g <- intersect(programs[[program]],rownames(z)); if(length(g)<3) next
    score <- colMeans(t(scale(t(z[g,,drop=FALSE]))),na.rm=TRUE)
    target <- if(program=="endothelial") "Endothelial cells" else "T cells & NK cells"
    idx <- !is.na(labels); score <- score[idx]; lab <- labels[idx]
    grp <- lab==target
    wt <- if(sum(grp)>=3 && sum(!grp)>=3) suppressWarnings(wilcox.test(score[grp],score[!grp],exact=FALSE)) else NULL
    med_target <- unname(median(score[grp],na.rm=TRUE)); med_other <- unname(median(score[!grp],na.rm=TRUE))
    rows <- append(rows,list(data.frame(dataset=unname("GSE111672"),sample=unname(sample_name),program=unname(program),target_label=unname(target),n_genes_used=length(g),n_target=unname(sum(grp)),n_other=unname(sum(!grp)),median_target=med_target,median_other=med_other,delta=unname(med_target-med_other),wilcox_P=if(is.null(wt)) NA_real_ else unname(wt$p.value),stringsAsFactors=FALSE,row.names=NULL)))
    avail <- append(avail,list(data.frame(dataset="GSE111672",sample=sample_name,program=program,n_genes_total=length(programs[[program]]),n_genes_found=length(g),genes_found=paste(g,collapse=";"),stringsAsFactors=FALSE)))
  }
}
if(length(rows)) fwrite(rbindlist(rows),file.path(out_dir,"GSE111672_program_validation.tsv"),sep="\t",na="NA")
if(length(avail)) fwrite(rbindlist(avail),file.path(out_dir,"GSE111672_program_gene_availability.tsv"),sep="\t",na="NA")
