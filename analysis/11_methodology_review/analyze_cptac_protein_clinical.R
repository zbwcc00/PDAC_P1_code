options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(survival)})
root <- "D:/PDAC_P1"; indir <- file.path(root,"data/05_protein/CPTAC_PAAD_cBioPortal"); outdir <- file.path(root,"analysis/11_methodology_review/CPTAC_PAAD")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
p <- fread(file.path(indir,"protein_quantification.tsv")); c <- fread(file.path(indir,"patient_clinical_data.tsv"))
map <- data.table(entrezGeneId=c(2034,6867,4082,9709,5175,7450,51705),gene=c("EPAS1","TACC1","MARCKS","HERPUD1","PECAM1","VWF","EMCN"))
p <- merge(p,map,by="entrezGeneId",all.x=TRUE)
prot <- p[,.(n_samples=.N,n_nonmissing=sum(!is.na(value)),median=median(value,na.rm=TRUE),IQR=IQR(value,na.rm=TRUE),mean=mean(value,na.rm=TRUE)),by=gene]
fwrite(prot,file.path(outdir,"CPTAC_protein_detection_summary.tsv"),sep="\t")
surv <- c[clinicalAttributeId %in% c("FOLLOW_UP_DAYS","VITAL_STATUS"),.(patientId,clinicalAttributeId,value)]
surv <- dcast(surv,patientId~clinicalAttributeId,value.var="value"); surv[,days:=as.numeric(FOLLOW_UP_DAYS)]; surv[,event:=as.integer(tolower(VITAL_STATUS)=="deceased")]; surv <- surv[is.finite(days)&days>0&!is.na(event)]
pm <- dcast(p[,.(patientId,gene,value)],patientId~gene,value.var="value",fun.aggregate=median)
dat <- merge(surv,pm,by="patientId"); genes <- intersect(c("EPAS1","TACC1","MARCKS","HERPUD1","PECAM1","VWF","EMCN"),names(dat)); rows <- list()
for(g in genes){if(sum(!is.na(dat[[g]]))<20) next; z<-dat[!is.na(get(g))]; med<-median(z[[g]],na.rm=TRUE); z[,high:=get(g)>=med]; fit<-coxph(Surv(days,event)~get(g),data=z); fit2<-coxph(Surv(days,event)~high,data=z); s<-summary(fit); s2<-summary(fit2); rows[[g]]<-data.frame(gene=g,n=nrow(z),events=sum(z$event),continuous_HR=exp(coef(fit)[1]),continuous_P=coef(s)[1,5],median_split_HR=exp(coef(fit2)[1]),median_split_P=coef(s2)[1,5],stringsAsFactors=FALSE)}
cox <- rbindlist(rows,fill=TRUE); fwrite(cox,file.path(outdir,"CPTAC_protein_survival_cox.tsv"),sep="\t",na="NA")
writeLines(c("# CPTAC-PAAD protein and clinical validation","",paste0("- cBioPortal study: paad_cptac_2021; 140 samples."),"- Protein profile: CPTAC protein quantification (log2-value); clinical endpoint uses follow-up days and vital status.",paste0("- Available protein records: ",nrow(p),"."),"- EPAS1 protein records are absent in this profile; HERPUD1 and the other marker proteins have partial or complete coverage.","- Survival associations are exploratory and unadjusted; they do not replace an independent, pre-specified prognostic validation cohort."),file.path(outdir,"CPTAC_protein_clinical_report.md"))
