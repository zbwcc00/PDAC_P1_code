options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1"
cand <- fread(file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
cand <- cand[major_lineage == "endothelial" & stable_high_confidence == TRUE]
g24 <- fread(file.path(root, "analysis/08_program_validation/GSE62452/GSE62452_paired_candidate_gene_results.tsv"))
g71 <- fread(file.path(root, "analysis/08_program_validation/GSE71729/GSE71729_primary_vs_normal_candidate_gene_results.tsv"))
setnames(g24, c("logFC", "P.Value", "adj.P.Val"), c("gse62452_logFC", "gse62452_P", "gse62452_FDR"))
setnames(g71, c("logFC", "P.Value", "adj.P.Val"), c("gse71729_logFC", "gse71729_P", "gse71729_FDR"))
res <- merge(cand[, .(gene_id, direction, median_logFC, mean_abs_logFC, n_FDR05)], g24[, .(gene, gse62452_logFC, gse62452_P, gse62452_FDR)], by.x="gene_id", by.y="gene", all.x=TRUE)
res <- merge(res, g71[, .(gene, gse71729_logFC, gse71729_P, gse71729_FDR)], by.x="gene_id", by.y="gene", all.x=TRUE)
res[, expected_sign := ifelse(direction == "tumor_up", 1, -1)]
res[, sign_62452 := sign(gse62452_logFC) == expected_sign]
res[, sign_71729 := sign(gse71729_logFC) == expected_sign]
res[, external_FDR_both := !is.na(gse62452_FDR) & !is.na(gse71729_FDR) & gse62452_FDR < 0.05 & gse71729_FDR < 0.05]
res[, external_direction_both := sign_62452 & sign_71729]
res[, target_priority := fifelse(external_FDR_both & external_direction_both, "Tier1", fifelse(external_direction_both, "Tier2", "Tier3"))]
res[, evidence_score := -log10(pmax(gse62452_FDR, 1e-300)) - log10(pmax(gse71729_FDR, 1e-300)) + abs(median_logFC)]
res[, priority_order := match(target_priority, c("Tier1", "Tier2", "Tier3"))]
res[is.na(evidence_score), evidence_score := -Inf]
setorder(res, priority_order, -evidence_score)
fwrite(res, file.path(root, "analysis/08_program_validation/endothelial_target_prioritization.tsv"), sep="\t")
sink(file.path(root, "analysis/08_program_validation/endothelial_target_prioritization_report.md")); cat("# Endothelial target prioritization\\n\\n"); cat("Tier 1 requires concordant direction and FDR<0.05 in both independent bulk cohorts (GSE62452 paired and GSE71729 primary-vs-normal).\\n\\n"); print(res); sink()
