options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1"
indir <- file.path(root,"data/03_functional/PRISM20Q2/raw")
outdir <- file.path(root,"analysis/11_methodology_review/PRISM20Q2")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
lfc <- fread(file.path(indir,"prism-repurposing-20q2-primary-screen-replicate-collapsed-logfold-change.csv"),check.names=FALSE)
setnames(lfc,1,"depmap_id")
cells <- fread(file.path(indir,"prism-repurposing-20q2-primary-screen-cell-line-info.csv"))
pdac <- cells[grepl("pancreas",primary_tissue,ignore.case=TRUE) | grepl("pancreas",secondary_tissue,ignore.case=TRUE),]
pdac_ids <- intersect(pdac$depmap_id, lfc$depmap_id)
fwrite(pdac,file.path(outdir,"PRISM20Q2_pancreas_cell_lines.tsv"),sep="\t")
mat <- as.matrix(lfc[match(pdac_ids,lfc$depmap_id),-1,with=FALSE]); rownames(mat) <- pdac_ids
treat <- fread(file.path(indir,"prism-repurposing-20q2-primary-screen-replicate-collapsed-treatment-info.csv"))
colkey <- colnames(mat)
# Match by the complete column_name to preserve the expression-matrix column order.
drug <- treat[match(colkey, column_name), .(broad_id, dose, screen_id, name, moa, target, disease.area, indication, smiles)]
drug[, column := colkey]
rows <- lapply(seq_len(ncol(mat)),function(j){v=mat[,j]; data.frame(column=colkey[j],broad_id=drug$broad_id[j],name=drug$name[j],moa=drug$moa[j],target=drug$target[j],screen_id=drug$screen_id[j],n_panc=sum(!is.na(v)),median_lfc=median(v,na.rm=TRUE),mean_lfc=mean(v,na.rm=TRUE),fraction_lfc_le_neg05=mean(v<=-0.5,na.rm=TRUE),stringsAsFactors=FALSE)})
summary <- rbindlist(rows,fill=TRUE)
summary <- summary[order(median_lfc)]
fwrite(summary,file.path(outdir,"PRISM20Q2_pancreas_drug_summary.tsv"),sep="\t",na="NA")
sel <- summary[grepl("Y-39983|triclabendazole",name,ignore.case=TRUE)]
fwrite(sel,file.path(outdir,"PRISM20Q2_EPAS1_candidate_drug_results.tsv"),sep="\t",na="NA")
top <- summary[n_panc>=20][1:min(50,.N)]
fwrite(top,file.path(outdir,"PRISM20Q2_top_pancreas_sensitive_drugs.tsv"),sep="\t",na="NA")
writeLines(c("# PRISM 20Q2 PDAC pharmacologic audit","",paste0("- Matrix dimensions: ",nrow(lfc)," cell lines × ",ncol(lfc)-1," treatment columns."),paste0("- Pancreas-labelled models: ",length(pdac_ids),"."),"- Effect metric: replicate-collapsed log-fold-change; more negative values indicate greater viability loss/sensitivity.",paste0("- Candidate-drug rows matching Y-39983 or triclabendazole: ",nrow(sel),"."),"- PRISM is an in-vitro orthogonal sensitivity layer; it does not establish endothelial specificity, clinical efficacy or target engagement."),file.path(outdir,"PRISM20Q2_pdac_report.md"))
