options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

root <- "D:/PDAC_P1"
raw <- file.path(root, "data/03_functional/DepMap24Q4/raw")
prism <- file.path(root, "data/03_functional/PRISM20Q2/raw")
out <- file.path(root, "analysis/11_methodology_review/DepMap_PRISM_pharmacology")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

expression_candidates <- c(file.path(raw, "OmicsExpressionProteinCodingGenesTPMLogp1.csv"), file.path(raw, "OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv"))
expr.file <- expression_candidates[file.exists(expression_candidates)][1]
if (is.na(expr.file)) stop("No DepMap protein-coding TPM matrix was found in: ", raw)
expression_type <- if (basename(expr.file) == "OmicsExpressionProteinCodingGenesTPMLogp1.csv") "non_stranded" else "stranded"

models <- fread(file.path(raw, "Model.csv"))
pdac <- models[OncotreeLineage == "Pancreas", .(ModelID, CellLineName, CCLEName, OncotreePrimaryDisease)]
candidates <- fread(file.path(root, "analysis/07_candidate_screen/candidate_genes_stable_high_confidence_all.tsv"))
up <- candidates[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_up", gene_id]
down <- candidates[major_lineage == "endothelial" & stable_high_confidence == TRUE & direction == "tumor_down", gene_id]
targets <- unique(c("EPAS1", "MARCKS", "HERPUD1", "TACC1", up, down, "PECAM1", "VWF", "EMCN"))

header <- names(fread(expr.file, nrows = 0, check.names = FALSE))
find_gene <- function(gene) { match <- header[startsWith(header, paste0(gene, " ("))]; if (length(match)) match[1] else NA_character_ }
expression_columns <- unique(na.omit(c(header[1], vapply(targets, find_gene, character(1)))))
expression <- fread(expr.file, select = expression_columns, check.names = FALSE)
setnames(expression, 1, "ModelID")
setnames(expression, names(expression), sub(" \\([0-9]+\\)$", "", names(expression)))

logfc.file <- file.path(prism, "prism-repurposing-20q2-primary-screen-replicate-collapsed-logfold-change.csv")
prism_header <- names(fread(logfc.file, nrows = 0, check.names = FALSE))
drug_columns <- c("BRD-K81916719-001-13-9::2.5::HTS", "BRD-K56751279-300-01-4::2.5::HTS")
stopifnot(all(drug_columns %in% prism_header))
drug_response <- fread(logfc.file, select = c(prism_header[1], drug_columns), check.names = FALSE)
setnames(drug_response, 1, "ModelID")
setnames(drug_response, drug_columns, c("triclabendazole_logFC", "Y39983_logFC"))

matched <- Reduce(function(left, right) merge(left, right, by = "ModelID"), list(pdac, expression, drug_response))
z_score <- function(values) as.numeric(scale(values))
vascular_genes <- intersect(c("PECAM1", "VWF", "EMCN"), names(matched))
if (length(vascular_genes)) matched[, vascular_marker_score := rowMeans(sapply(vascular_genes, function(gene) z_score(get(gene))), na.rm = TRUE)]
up_genes <- intersect(up, names(matched)); down_genes <- intersect(down, names(matched))
if (length(up_genes) && length(down_genes)) matched[, endothelial_program := rowMeans(sapply(up_genes, function(gene) z_score(get(gene))), na.rm = TRUE) - rowMeans(sapply(down_genes, function(gene) z_score(get(gene))), na.rm = TRUE)]

features <- intersect(c("EPAS1", "MARCKS", "HERPUD1", "TACC1", "endothelial_program", "vascular_marker_score"), names(matched))
results <- rbindlist(lapply(c("triclabendazole_logFC", "Y39983_logFC"), function(drug) rbindlist(lapply(features, function(feature) {
  subset <- matched[is.finite(get(feature)) & is.finite(get(drug))]
  if (nrow(subset) < 20) return(data.table(drug = drug, feature = feature, n = nrow(subset), spearman_rho = NA_real_, P = NA_real_, analysis_status = "below_prespecified_n20_gate"))
  test <- cor.test(subset[[feature]], subset[[drug]], method = "spearman", exact = FALSE)
  data.table(drug = drug, feature = feature, n = nrow(subset), spearman_rho = unname(test$estimate), P = test$p.value, analysis_status = "eligible_association_test")
}))))
results[, FDR := p.adjust(P, "BH")]

fwrite(matched, file.path(out, "DepMap_PRISM_PDAC_matched_features.tsv"), sep = "\t")
fwrite(results, file.path(out, "DepMap_PRISM_association_results.tsv"), sep = "\t")
writeLines(c("# DepMap expression-PRISM pharmacologic association audit", "", paste0("- Expression input: `", basename(expr.file), "` (", expression_type, ")."), paste0("- PDAC models in DepMap Model.csv: ", nrow(pdac), "; matched PRISM/expression models: ", nrow(matched), "."), "- Prespecified gate: a drug-feature association requires >=20 PDAC models with finite expression and logFC; otherwise no correlation, regression, or machine-learning conclusion is made.", "- PRISM logFC: more negative values indicate greater in-vitro viability loss.", "- This is hypothesis-generating tumor-cell pharmacogenomics, not direct target engagement, endothelial specificity, or clinical efficacy."), file.path(out, "DepMap_PRISM_pharmacology_report.md"))
