options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(scTenifoldKnk); library(data.table)})

set.seed(20260906)
root <- "D:/PDAC_P1"
outdir <- file.path(root, "analysis/17_P0_reinforcement")
resdir <- file.path(root, "analysis/13_virtual_perturbation/results")
objdir <- file.path(root, "analysis/03_scrna_annotation/seurat_annotated")
dir.create(file.path(outdir, "empirical_null_results"), showWarnings = FALSE)

datasets <- c(
  GSE154778_primary = "GSE154778__primary_tumor.rds",
  GSE155698_primary = "GSE155698__primary_tumor.rds",
  GSE212966_primary = "GSE212966__primary_tumor.rds"
)
downstream_programs <- list(
  endothelial_downstream = c("CLDN5", "KDR", "PECAM1", "VWF", "EMCN", "ESAM", "RAMP2", "PLVAP", "ENG"),
  angiogenesis_downstream = c("VEGFA", "KDR", "PECAM1", "VWF", "EMCN", "ENG")
)

summarize_response <- function(diff_regulation, dataset, target, target_class) {
  rbindlist(lapply(names(downstream_programs), function(program) {
    genes <- intersect(downstream_programs[[program]], diff_regulation$gene)
    d <- diff_regulation[match(genes, diff_regulation$gene), ]
    data.table(
      dataset = dataset,
      target = target,
      target_class = target_class,
      endpoint = program,
      n_endpoint_genes_tested = length(genes),
      n_FDR05 = sum(d$p.adj < 0.05, na.rm = TRUE),
      downstream_network_score = sum(-log10(pmax(d$p.adj, 1e-300)), na.rm = TRUE),
      downstream_mean_abs_Z = mean(abs(d$Z), na.rm = TRUE)
    )
  }))
}

all_matches <- list()
all_scores <- list()
for (dataset in names(datasets)) {
  message("Preparing ", dataset)
  epas1_rds <- readRDS(file.path(resdir, paste0(dataset, "__endothelial__EPAS1__KO.rds")))
  wt <- epas1_rds$tensorNetworks$WT
  network_genes <- rownames(wt)
  degree <- rowSums(abs(wt))
  
  object <- readRDS(file.path(objdir, datasets[[dataset]]))
  assay <- if ("RNA" %in% names(object@assays)) "RNA" else DefaultAssay(object)
  counts <- GetAssayData(object, assay = assay, layer = "counts")
  metadata <- object@meta.data
  endothelial_cells <- rownames(metadata)[metadata$major_lineage == "endothelial"]
  sub <- counts[intersect(network_genes, rownames(counts)), endothelial_cells, drop = FALSE]
  sub <- sub[, Matrix::colSums(sub) >= 100, drop = FALSE]
  sub <- sub[network_genes[network_genes %in% rownames(sub)], , drop = FALSE]
  log_sub <- log1p(t(t(sub) / pmax(Matrix::colSums(sub), 1)) * 1e4)
  
  properties <- data.table(
    gene = rownames(sub),
    detection = Matrix::rowMeans(sub > 0),
    mean_log_expression = Matrix::rowMeans(log_sub),
    network_degree = degree[rownames(sub)]
  )
  epas1_properties <- properties[gene == "EPAS1"]
  excluded <- unique(c("EPAS1", unlist(downstream_programs), "TACC1", "MARCKS", "HERPUD1"))
  controls <- properties[!gene %in% excluded & detection > 0.02 & is.finite(network_degree)]
  for (column in c("detection", "mean_log_expression", "network_degree")) {
    scale_value <- sd(controls[[column]], na.rm = TRUE)
    if (!is.finite(scale_value) || scale_value == 0) scale_value <- 1
    controls[[paste0("z_", column)]] <- (controls[[column]] - epas1_properties[[column]]) / scale_value
  }
  controls[, matching_distance := sqrt(z_detection^2 + z_mean_log_expression^2 + z_network_degree^2)]
  setorder(controls, matching_distance)
  matched <- controls[seq_len(min(20L, .N))]
  matched[, c("dataset", "n_endothelial_cells", "matching_rule") := list(
    dataset, ncol(sub),
    "Nearest 20 controls after exclusion of EPAS1, endpoint genes, and other prioritized targets; Euclidean distance on detection, mean log expression, and WT network degree."
  )]
  all_matches[[dataset]] <- matched
  
  epas1_score <- summarize_response(epas1_rds$diffRegulation, dataset, "EPAS1", "observed")
  all_scores[[length(all_scores) + 1L]] <- epas1_score
  
  for (control_gene in matched$gene) {
    message("  KO null: ", control_gene)
    result_path <- file.path(outdir, "empirical_null_results", paste0(dataset, "__endothelial__", control_gene, "__KO.rds"))
    control_result <- tryCatch({
      if (file.exists(result_path)) readRDS(result_path) else {
        value <- scTenifoldKnk(sub, qc = FALSE, gKO = control_gene, nc_nNet = 2,
          nc_nCells = min(100, ncol(sub)), nc_nComp = 3, td_K = 3, ma_nDim = 2, nCores = 1)
        saveRDS(value, result_path, compress = "xz")
        value
      }
    }, error = function(e) e)
    if (inherits(control_result, "error")) {
      all_scores[[length(all_scores) + 1L]] <- data.table(
        dataset = dataset, target = control_gene, target_class = "matched_null_error",
        endpoint = NA_character_, n_endpoint_genes_tested = NA_integer_, n_FDR05 = NA_integer_,
        downstream_network_score = NA_real_, downstream_mean_abs_Z = NA_real_,
        error_message = conditionMessage(control_result)
      )
    } else {
      all_scores[[length(all_scores) + 1L]] <- summarize_response(control_result$diffRegulation, dataset, control_gene, "matched_null")
      rm(control_result); gc()
    }
  }
  rm(object, counts, sub, log_sub, epas1_rds); gc()
}

matches <- rbindlist(all_matches, fill = TRUE)
scores <- rbindlist(all_scores, fill = TRUE)
fwrite(matches, file.path(outdir, "EPAS1_KO_matched_null_gene_properties.tsv"), sep = "\t")
fwrite(scores, file.path(outdir, "EPAS1_KO_empirical_null_program_responses.tsv"), sep = "\t")

empirical <- scores[target_class == "observed"][, {
  null_values <- scores[target_class == "matched_null" & dataset == .BY$dataset & endpoint == .BY$endpoint, downstream_network_score]
  null_values <- null_values[is.finite(null_values)]
  observed <- downstream_network_score
  data.table(
    observed_downstream_network_score = observed,
    n_matched_null = length(null_values),
    n_null_at_least_observed = sum(null_values >= observed),
    empirical_one_sided_p = (1 + sum(null_values >= observed)) / (1 + length(null_values)),
    observed_rank_among_nulls = 1 + sum(null_values > observed)
  )
}, by = .(dataset, endpoint)]
fwrite(empirical, file.path(outdir, "EPAS1_KO_empirical_null_summary.tsv"), sep = "\t")

readme <- c(
  "# EPAS1 scTenifoldKnk empirical-null analysis",
  "",
  "Endpoint programs exclude EPAS1 itself, so a target gene cannot create a trivial positive score.",
  "Each cohort uses 20 non-prioritized endothelial genes matched on detection rate, mean log expression, and WT-network weighted degree.",
  "The empirical p value is (1 + number of null scores >= observed score) / (1 + number of completed null scores).",
  "This benchmark evaluates network-model specificity only. It remains an in silico perturbation and does not establish genetic causality."
)
writeLines(readme, file.path(outdir, "EPAS1_KO_empirical_null_readme.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "03_sessionInfo.txt"))

