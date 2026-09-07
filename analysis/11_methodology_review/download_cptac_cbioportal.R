suppressPackageStartupMessages({library(jsonlite);library(httr)})
root <- "D:/PDAC_P1"
out <- file.path(root,"data/05_protein/CPTAC_PAAD_cBioPortal")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
getj <- function(u) fromJSON(u, simplifyDataFrame=TRUE)
study <- "paad_cptac_2021"
samples <- getj(paste0("https://www.cbioportal.org/api/studies/",study,"/samples?projection=SUMMARY"))
ids <- samples$sampleId
body <- list(entrezGeneIds=as.integer(c(2034,6867,4082,9709,5175,7450,51705)), sampleIds=as.character(ids))
resp <- httr::POST(paste0("https://www.cbioportal.org/api/molecular-profiles/",study,"_protein_quantification/molecular-data/fetch"), body=body, encode="json")
stop_for_status(resp)
protein <- content(resp,as="text",encoding="UTF-8") |> fromJSON()
write.table(protein,file.path(out,"protein_quantification.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
clinical <- getj(paste0("https://www.cbioportal.org/api/studies/",study,"/clinical-data?clinicalDataType=PATIENT&projection=DETAILED"))
write.table(clinical,file.path(out,"patient_clinical_data.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(samples,file.path(out,"sample_metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("samples",nrow(samples),"protein",length(protein),"clinical",nrow(clinical),"\n")
