options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1"
run_dir <- file.path(root,"analysis/04_cnv_copykat/GSE154778")
inp <- readRDS(file.path(run_dir,"copykat_input_sampled.rds"))
pred <- fread(file.path(run_dir,"PDAC_GSE154778_copykat_prediction.txt"))
setnames(pred,c("cell.names","copykat.pred"),c("cell","copykat_prediction"))
md <- as.data.table(inp$metadata)
calls <- merge(pred,md[,.(cell,dataset,gsm,tissue_context,patient_id,cnv_role)],by="cell",all.x=TRUE)
calls[, predicted_aneuploid := tolower(copykat_prediction)=="aneuploid"]
fwrite(calls,file.path(run_dir,"copykat_cell_calls.tsv"),sep="\t",na="NA")
summary <- calls[,.(sampled_cells=.N,aneuploid_cells=sum(predicted_aneuploid,na.rm=TRUE),aneuploid_fraction=mean(predicted_aneuploid,na.rm=TRUE)),by=.(dataset,gsm,tissue_context,patient_id,cnv_role)]
fwrite(summary,file.path(run_dir,"copykat_sample_summary.tsv"),sep="\t",na="NA")
fwrite(data.table(dataset="GSE154778",status="core_prediction_complete",sampled_cells=nrow(calls),aneuploid_fraction=mean(calls$predicted_aneuploid,na.rm=TRUE)),file.path(run_dir,"copykat_status.tsv"),sep="\t",na="NA")
