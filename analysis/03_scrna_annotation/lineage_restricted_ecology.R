options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1"
in_dir <- file.path(root, "analysis/03_scrna_annotation")
cells <- fread(file.path(in_dir, "PDAC_scRNA_cell_major_lineage_scores.tsv"))
flags <- fread(file.path(root, "analysis/02_scrna_qc/PDAC_scRNA_sample_review_flags.tsv"))
cells[, is_malignant := major_lineage == "malignant_epithelial"]
cells[, is_caf := major_lineage == "fibroblast_CAF"]
cells[, is_myeloid := major_lineage == "myeloid"]
cells[, is_cd8 := major_lineage == "CD8"]
sample <- cells[, .(
  total_singlet_cells = .N,
  malignant_prop = mean(is_malignant), CAF_prop = mean(is_caf), myeloid_prop = mean(is_myeloid), CD8_prop = mean(is_cd8),
  malignant_score_all = mean(module_malignant_epithelial_program), CAF_score_all = mean(module_CAF_program),
  myeloid_score_all = mean(module_myeloid_program), CD8_activation_all = mean(module_CD8_activation), CD8_exhaustion_all = mean(module_CD8_exhaustion),
  malignant_score_restricted = ifelse(sum(is_malignant)>0, mean(module_malignant_epithelial_program[is_malignant]), NA_real_),
  CAF_score_restricted = ifelse(sum(is_caf)>0, mean(module_CAF_program[is_caf]), NA_real_),
  myeloid_score_restricted = ifelse(sum(is_myeloid)>0, mean(module_myeloid_program[is_myeloid]), NA_real_),
  CD8_activation_restricted = ifelse(sum(is_cd8)>0, mean(module_CD8_activation[is_cd8]), NA_real_),
  CD8_exhaustion_restricted = ifelse(sum(is_cd8)>0, mean(module_CD8_exhaustion[is_cd8]), NA_real_)
), by = .(dataset, gsm, tissue_context, patient_id)]
sample <- merge(sample, flags[, .(dataset, gsm, review_flag)], by=c("dataset","gsm"), all.x=TRUE)
sample[is.na(review_flag), review_flag := "acceptable"]
fwrite(sample, file.path(in_dir, "PDAC_scRNA_patient_lineage_restricted_ecology.tsv"), sep="\t", na="NA")

safe_cor <- function(a,b) if(sum(complete.cases(a,b))>=5 && length(unique(a[complete.cases(a,b)]))>2 && length(unique(b[complete.cases(a,b)]))>2) suppressWarnings(cor(a,b,method="spearman",use="pairwise.complete.obs")) else NA_real_
run_cor <- function(dat, analysis_set) {
  out <- list()
  for(ds in unique(dat$dataset)) {
    x <- dat[dataset==ds & tissue_context=="primary_tumor"]
    if(nrow(x)<6) next
    out <- append(out,list(data.frame(analysis_set=analysis_set,dataset=ds,n_patients=nrow(x),
      rho_malignant_restricted_CAF_prop=safe_cor(x$malignant_score_restricted,x$CAF_prop),
      rho_malignant_restricted_myeloid_prop=safe_cor(x$malignant_score_restricted,x$myeloid_prop),
      rho_malignant_restricted_CD8_exhaustion=safe_cor(x$malignant_score_restricted,x$CD8_exhaustion_restricted),
      rho_CAF_prop_CD8_prop=safe_cor(x$CAF_prop,x$CD8_prop),
      rho_myeloid_prop_CD8_prop=safe_cor(x$myeloid_prop,x$CD8_prop), stringsAsFactors=FALSE)))
  }
  if(length(out)) rbindlist(out) else data.table()
}
cor_all <- run_cor(sample, "all_primary")
cor_sens <- run_cor(sample[review_flag != "exclude_or_manual_review"], "exclude_extreme_QC")
fwrite(rbind(cor_all, cor_sens), file.path(in_dir, "PDAC_scRNA_lineage_restricted_correlations.tsv"), sep="\t", na="NA")
