options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1"
pr <- fread(file.path(root,"analysis/11_methodology_review/PRISM20Q2/PRISM20Q2_EPAS1_candidate_drug_results.tsv"))
st <- fread(file.path(root,"analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/epas1_final_structure_evidence.tsv"))
out <- merge(st,pr[,.(cmap_name=name,prism_n_panc=n_panc,prism_median_lfc=median_lfc,prism_mean_lfc=mean_lfc,prism_fraction_lfc_le_neg05=fraction_lfc_le_neg05)],by="cmap_name",all.x=TRUE)
out[,integrated_decision:=fifelse(cmap_name=="triclabendazole","retain_primary_prism_supported", "retain_as_negative_pharmacologic_comparator")]
fwrite(out,file.path(root,"analysis/11_methodology_review/PRISM20Q2/EPAS1_structure_prism_integrated.tsv"),sep="\t",na="NA")
writeLines(c("# EPAS1 structure–PRISM integrated decision","", "- Triclabendazole: DrugReflector recurrent, five-conformer docking stable, redocking QC passed, and PRISM median logFC −0.856 across 33 pancreas models; retain as the primary translational candidate.", "- Y-39983: DrugReflector recurrent and structurally stable, but PRISM median logFC −0.095 across 36 pancreas models with only 16.7% below −0.5; retain only as a negative pharmacologic comparator and do not claim PDAC sensitivity.", "- PRISM viability does not establish EPAS1 target engagement or endothelial specificity."),file.path(root,"analysis/11_methodology_review/PRISM20Q2/EPAS1_structure_prism_integrated_report.md"))
