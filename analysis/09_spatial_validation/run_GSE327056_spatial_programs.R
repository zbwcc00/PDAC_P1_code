options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(Matrix);library(data.table);library(edgeR);library(ggplot2)})
root <- "D:/PDAC_P1"
raw <- file.path(root,"data/03_spatial/GSE327056/raw")
out <- file.path(root,"analysis/09_spatial_validation")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"spot_tables"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"plots"),recursive=TRUE,showWarnings=FALSE)
cand <- fread(file.path(root,"analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))[candidate_status=="screening_candidate"]
programs <- split(cand$gene_id,cand$major_lineage)
programs <- programs[names(programs) %in% c("endothelial","T_NK")]
samples <- data.table(tag=c("A1","B1","C1","D1"),gsm=c("GSM9647219","GSM9647220","GSM9647221","GSM9647222"),context=c("adjacent_tumor","tumor","tumor_stroma","normal_pancreas"))
all_spots <- list(); sample_sum <- list(); gene_avail <- list()
for(i in seq_len(nrow(samples))) {
  s <- samples[i]; mfile <- file.path(raw,paste0(s$gsm,"_",s$tag,"_matrix.mtx.gz")); ffile <- file.path(raw,paste0(s$gsm,"_",s$tag,"_features.tsv.gz")); bfile <- file.path(raw,paste0(s$gsm,"_",s$tag,"_barcodes.tsv.gz")); pfile <- file.path(raw,paste0(s$gsm,"_",s$tag,"_tissue_positions.csv.gz"))
  mat <- readMM(gzfile(mfile)); feat <- fread(ffile,header=FALSE); bar <- fread(bfile,header=FALSE); pos <- fread(pfile)
  symbols <- if(ncol(feat)>=2) feat[[2]] else feat[[1]]; rownames(mat) <- make.unique(symbols); colnames(mat) <- bar[[1]]
  pos <- pos[barcode %in% colnames(mat) & in_tissue==1]; ids <- intersect(pos$barcode,colnames(mat)); pos <- pos[match(ids,pos$barcode)]; mat <- mat[,ids,drop=FALSE]
  y <- DGEList(counts=mat); logcpm <- cpm(y,log=TRUE,prior.count=2)
  one <- data.table(dataset="GSE327056",tag=s$tag,gsm=s$gsm,context=s$context,barcode=ids,array_row=pos$array_row,array_col=pos$array_col,pxl_row_in_fullres=pos$pxl_row_in_fullres,pxl_col_in_fullres=pos$pxl_col_in_fullres)
  for(program in names(programs)) {
    g <- intersect(programs[[program]],rownames(logcpm)); score <- if(length(g)>=3) colMeans(t(scale(t(logcpm[g,,drop=FALSE]))),na.rm=TRUE) else rep(NA_real_,ncol(logcpm)); one[[paste0(program,"_score")]] <- score; gene_avail <- append(gene_avail,list(data.frame(tag=s$tag,context=s$context,program=program,n_genes_total=length(programs[[program]]),n_genes_found=length(g),genes_found=paste(g,collapse=";"),stringsAsFactors=FALSE)))
    sample_sum <- append(sample_sum,list(data.frame(tag=s$tag,gsm=s$gsm,context=s$context,program=program,n_spots=ncol(logcpm),median_score=median(score,na.rm=TRUE),mean_score=mean(score,na.rm=TRUE),q25=quantile(score,0.25,na.rm=TRUE),q75=quantile(score,0.75,na.rm=TRUE),stringsAsFactors=FALSE)))
    plot_dt <- one[is.finite(get(paste0(program,"_score")))]; plot_dt[,score_value:=get(paste0(program,"_score"))]
    p <- ggplot(plot_dt,aes(x=pxl_col_in_fullres,y=pxl_row_in_fullres,color=score_value)) + geom_point(size=0.55) + scale_y_reverse() + coord_fixed() + scale_color_viridis_c(option="magma") + labs(title=paste0("GSE327056 ",s$tag," ",s$context," ",program),x="pixel column",y="pixel row",color="program score") + theme_void() + theme(plot.title=element_text(size=10))
    ggsave(file.path(out,"plots",paste0(s$tag,"_",program,"_spatial_score.png")),p,width=6,height=5,dpi=220)
  }
  fwrite(one,file.path(out,"spot_tables",paste0(s$tag,"_spot_program_scores.tsv")),sep="\t",na="NA"); all_spots <- append(all_spots,list(one)); rm(mat,logcpm,one); gc(verbose=FALSE)
}
fwrite(rbindlist(all_spots,fill=TRUE),file.path(out,"GSE327056_all_spot_program_scores.tsv"),sep="\t",na="NA")
fwrite(rbindlist(sample_sum,fill=TRUE),file.path(out,"GSE327056_spatial_program_summary.tsv"),sep="\t",na="NA")
fwrite(rbindlist(gene_avail,fill=TRUE),file.path(out,"GSE327056_spatial_gene_availability.tsv"),sep="\t",na="NA")
writeLines(capture.output(sessionInfo()),file.path(out,"spatial_program_sessionInfo.txt"))
