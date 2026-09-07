options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1/analysis/03_scrna_annotation"
cells <- fread(file.path(root,"PDAC_scRNA_cell_major_lineage_scores.tsv"))
cells[dataset=="GSE155698", patient_id := sub("[A-Z]$","",patient_id)]
fwrite(cells,file.path(root,"PDAC_scRNA_cell_major_lineage_scores.tsv"),sep="\t",na="NA")
denom <- cells[,.(total_singlet_cells=.N),by=.(dataset,gsm,tissue_context,patient_id)]
counts <- dcast(cells,dataset+gsm+tissue_context+patient_id~major_lineage,fun.aggregate=length,value.var="cell",fill=0)
wide <- merge(denom,counts,by=c("dataset","gsm","tissue_context","patient_id"),all=TRUE)
lineage_cols <- setdiff(names(wide),c("dataset","gsm","tissue_context","patient_id","total_singlet_cells"))
for(nm in lineage_cols) wide[[paste0("prop_",nm)]] <- wide[[nm]]/wide$total_singlet_cells
mods <- cells[,.(malignant_epithelial_program=mean(module_malignant_epithelial_program),CAF_program=mean(module_CAF_program),myeloid_program=mean(module_myeloid_program),CD8_activation=mean(module_CD8_activation),CD8_exhaustion=mean(module_CD8_exhaustion)),by=.(dataset,gsm,tissue_context,patient_id)]
wide <- merge(wide,mods,by=c("dataset","gsm","tissue_context","patient_id"),all=TRUE)
fwrite(wide,file.path(root,"PDAC_scRNA_patient_lineage_summary.tsv"),sep="\t",na="NA")
