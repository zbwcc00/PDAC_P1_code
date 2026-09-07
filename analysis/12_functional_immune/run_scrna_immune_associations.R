options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1"
out_dir <- file.path(root,"analysis/12_functional_immune")
sc <- fread(file.path(root,"analysis/08_program_validation/scRNA_patient_program_scores.tsv"))
cell <- fread(file.path(root,"analysis/03_scrna_annotation/PDAC_scRNA_cell_major_lineage_scores.tsv"))
endo <- sc[program=="endothelial" & major_lineage=="endothelial" & analysis_eligible==TRUE]
endo <- endo[,.(endothelial_score=median(program_score,na.rm=TRUE),n_endothelial_groups=.N),by=.(dataset,patient_id,tissue_context)]
other <- cell[,.(CAF_score=median(module_CAF_program,na.rm=TRUE),myeloid_score=median(module_myeloid_program,na.rm=TRUE),CD8_activation=median(module_CD8_activation,na.rm=TRUE),CD8_exhaustion=median(module_CD8_exhaustion,na.rm=TRUE),n_cells=.N),by=.(dataset,patient_id,tissue_context)]
merged <- merge(endo,other,by=c("dataset","patient_id","tissue_context"),all=FALSE)
fwrite(merged,file.path(out_dir,"scRNA_patient_endothelial_immune_merged.tsv"),sep="\t")
vars <- c("CAF_score","myeloid_score","CD8_activation","CD8_exhaustion")
rows <- lapply(vars,function(v){ok <- complete.cases(merged$endothelial_score,merged[[v]]); if(sum(ok)<4) return(data.frame(variable=v,n=sum(ok),rho=NA_real_,p=NA_real_)); ct <- suppressWarnings(cor.test(merged$endothelial_score[ok],merged[[v]][ok],method="spearman",exact=FALSE)); data.frame(variable=v,n=sum(ok),rho=unname(ct$estimate),p=ct$p.value)})
fwrite(rbindlist(rows),file.path(out_dir,"scRNA_endothelial_immune_correlations.tsv"),sep="\t")
by_ds <- merged[,lapply(vars,function(v){ok<-complete.cases(endothelial_score,get(v));if(sum(ok)<4)return(NA_real_);unname(cor(endothelial_score[ok],get(v)[ok],method="spearman"))}),by=.(dataset,tissue_context)]
setnames(by_ds,old=names(by_ds)[-(1:2)],new=paste0(vars,"_rho"))
fwrite(by_ds,file.path(out_dir,"scRNA_endothelial_immune_correlations_by_dataset.tsv"),sep="\t")
