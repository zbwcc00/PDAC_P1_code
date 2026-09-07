options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(Matrix);library(edgeR)})
root <- "D:/PDAC_P1"
indir <- file.path(root,"data/02_scrna_external/GSE335452/raw/spatial")
outdir <- file.path(root,"analysis/09_spatial_validation/GSE335452")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
cand <- fread(file.path(root,"analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
endo_up <- cand[major_lineage=="endothelial" & stable_high_confidence==TRUE & direction=="tumor_up",gene_id]
endo_down <- cand[major_lineage=="endothelial" & stable_high_confidence==TRUE & direction=="tumor_down",gene_id]
score_signed <- function(z,up,down) {u<-intersect(up,rownames(z));d<-intersect(down,rownames(z));a<-if(length(u))colMeans(z[u,,drop=FALSE])else rep(0,ncol(z));b<-if(length(d))colMeans(z[d,,drop=FALSE])else rep(0,ncol(z)); (a-b)/(as.integer(length(u)>0)+as.integer(length(d)>0))}
prefix <- sub("[.]matrix[.]mtx[.]gz$","",basename(list.files(indir,pattern="matrix.mtx.gz")))
rows <- list(); sample_rows <- list()
for (pr in prefix) {
  feat <- fread(file.path(indir,paste0(pr,".features.tsv.gz")),header=FALSE); genes<-as.character(feat[[2]])
  mat <- readMM(gzfile(file.path(indir,paste0(pr,".matrix.mtx.gz")))); rownames(mat)<-make.unique(genes)
  lcp <- cpm(DGEList(counts=mat),log=TRUE,prior.count=2)
  # leave EPAS1 out of the endothelial signature before co-localization testing.
  sc <- score_signed(t(scale(t(lcp))),setdiff(endo_up,"EPAS1"),setdiff(endo_down,"EPAS1"))
  ep <- lcp["EPAS1",]; vm <- colMeans(t(scale(t(lcp[intersect(c("PECAM1","VWF","EMCN"),rownames(lcp)),,drop=FALSE]))),na.rm=TRUE)
  q <- quantile(ep,c(.25,.75),na.rm=TRUE); high<-ep>=q[2]; low<-ep<=q[1]
  one <- data.frame(sample=pr,n_spots=ncol(lcp),EPAS1_program_rho=cor(ep,sc,method="spearman"),EPAS1_marker_rho=cor(ep,vm,method="spearman"),program_marker_rho=cor(sc,vm,method="spearman"),program_delta_high_vs_low=median(sc[high])-median(sc[low]),program_wilcox_P=wilcox.test(sc[high],sc[low],exact=FALSE)$p.value,stringsAsFactors=FALSE)
  rows[[pr]] <- one
  sample_rows[[pr]] <- data.frame(sample=pr,spot=seq_len(ncol(lcp)),EPAS1_logCPM=ep,endothelial_program_excluding_EPAS1=sc,vascular_marker_score=vm,EPAS1_quartile=ifelse(high,"high",ifelse(low,"low","middle")))
}
res <- rbindlist(rows); res[,`:=`(rho_meta_mean=mean(EPAS1_program_rho),rho_meta_median=median(EPAS1_program_rho))]
fwrite(res,file.path(outdir,"GSE335452_patient_spatial_anchor_results.tsv"),sep="\t")
fwrite(rbindlist(sample_rows),file.path(outdir,"GSE335452_spot_anchor_scores.tsv"),sep="\t")
writeLines(c("# GSE335452 multi-patient spatial anchor validation","", "- Five PDAC Visium-like matrices (w1v, x2v, c3v, l4v, z5v) were analyzed independently.", "- The signed endothelial program was calculated after removing EPAS1; thus the EPAS1-program correlation is not a self-correlation.", "- EPAS1-high and -low spots use within-patient upper/lower quartiles.", "- GEO supplementary data did not include tissue-position/image files, so this analysis supports molecular co-localization across spatial spots but cannot produce histology-aligned maps.", "- Spot-level P values are descriptive; patient-level replication across all five matrices is the principal evidence."),file.path(outdir,"GSE335452_spatial_anchor_report.md"))
