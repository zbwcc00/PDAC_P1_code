options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

root <- "D:/PDAC_P1"
input_file <- file.path(root, "analysis/11_methodology_review/DepMap_PRISM_pharmacology/DepMap_PRISM_PDAC_matched_features.tsv")
out <- file.path(root, "analysis/14_prediction/DepMap_PRISM_drug_ML")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

set.seed(20260903)
data <- fread(input_file)
features <- intersect(c("EPAS1", "MARCKS", "HERPUD1", "TACC1", "endothelial_program", "vascular_marker_score"), names(data))
drugs <- intersect(c("triclabendazole_logFC", "Y39983_logFC"), names(data))

cross_validate_drug <- function(drug, repetitions = 100, folds = 5) {
  complete <- data[complete.cases(data[, c(features, drug), with = FALSE])]
  if (nrow(complete) < 20) return(list(metrics = data.table(), selections = data.table(), n = nrow(complete)))
  outcome <- complete[[drug]]
  predictor_matrix <- as.matrix(complete[, ..features])
  metrics <- vector("list", repetitions)
  selections <- vector("list", repetitions * folds)
  selection_index <- 1L

  for (repetition in seq_len(repetitions)) {
    fold_id <- sample(rep(seq_len(folds), length.out = nrow(complete)))
    predicted <- rep(NA_real_, nrow(complete))
    baseline <- rep(NA_real_, nrow(complete))
    for (fold in seq_len(folds)) {
      test_index <- which(fold_id == fold)
      train_index <- which(fold_id != fold)
      train_x <- predictor_matrix[train_index, , drop = FALSE]
      test_x <- predictor_matrix[test_index, , drop = FALSE]
      center <- colMeans(train_x)
      scale_value <- apply(train_x, 2, sd)
      scale_value[scale_value == 0 | !is.finite(scale_value)] <- 1
      train_x <- sweep(sweep(train_x, 2, center, "-"), 2, scale_value, "/")
      test_x <- sweep(sweep(test_x, 2, center, "-"), 2, scale_value, "/")
      inner_folds <- min(5, length(train_index))
      fitted <- cv.glmnet(train_x, outcome[train_index], alpha = 0.5, nfolds = inner_folds, standardize = FALSE, family = "gaussian")
      predicted[test_index] <- as.numeric(predict(fitted, newx = test_x, s = "lambda.1se"))
      baseline[test_index] <- mean(outcome[train_index])
      coefficients <- as.matrix(coef(fitted, s = "lambda.1se"))[-1, 1]
      selections[[selection_index]] <- data.table(drug = drug, repetition = repetition, fold = fold, feature = names(coefficients), selected = coefficients != 0, coefficient = coefficients)
      selection_index <- selection_index + 1L
    }
    total_ss <- sum((outcome - mean(outcome))^2)
    metrics[[repetition]] <- data.table(
      drug = drug,
      repetition = repetition,
      n = nrow(complete),
      elastic_net_R2 = 1 - sum((outcome - predicted)^2) / total_ss,
      baseline_R2 = 1 - sum((outcome - baseline)^2) / total_ss,
      elastic_net_RMSE = sqrt(mean((outcome - predicted)^2)),
      baseline_RMSE = sqrt(mean((outcome - baseline)^2)),
      spearman_rho = unname(cor(outcome, predicted, method = "spearman"))
    )
  }
  list(metrics = rbindlist(metrics), selections = rbindlist(selections), n = nrow(complete))
}

results <- lapply(drugs, cross_validate_drug)
metric_table <- rbindlist(lapply(results, `[[`, "metrics"), fill = TRUE)
selection_table <- rbindlist(lapply(results, `[[`, "selections"), fill = TRUE)
summary_table <- metric_table[, .(
  n = unique(n),
  repetitions = .N,
  median_elastic_net_R2 = median(elastic_net_R2),
  R2_2.5pct = quantile(elastic_net_R2, 0.025),
  R2_97.5pct = quantile(elastic_net_R2, 0.975),
  median_baseline_R2 = median(baseline_R2),
  median_elastic_net_RMSE = median(elastic_net_RMSE),
  median_baseline_RMSE = median(baseline_RMSE),
  median_spearman_rho = median(spearman_rho),
  proportion_R2_above_zero = mean(elastic_net_R2 > 0)
), by = drug]
summary_table[, model_decision := ifelse(median_elastic_net_R2 > 0.05 & R2_2.5pct > 0, "exploratory_signal_only", "do_not_retain")]
selection_summary <- selection_table[, .(selection_frequency = mean(selected), median_nonzero_coefficient = median(coefficient[selected], na.rm = TRUE)), by = .(drug, feature)]

fwrite(metric_table, file.path(out, "drug_ml_repeated_CV_metrics.tsv"), sep = "\t")
fwrite(summary_table, file.path(out, "drug_ml_repeated_CV_summary.tsv"), sep = "\t")
fwrite(selection_summary, file.path(out, "drug_ml_feature_stability.tsv"), sep = "\t")
writeLines(c(
  "# DepMap PRISM drug-response machine-learning sensitivity analysis",
  "",
  "- Model: Elastic Net Gaussian regression with alpha fixed at 0.5 and lambda.1se selected in the inner cross-validation loop.",
  "- Predictors were prespecified: EPAS1, MARCKS, HERPUD1, TACC1, endothelial program, and vascular-marker score.",
  "- Evaluation: 100 repetitions of five-fold outer cross-validation, compared with a training-fold mean-response baseline.",
  "- Retention gate: median outer-CV R2 > 0.05 and its 2.5th percentile > 0. This stringent gate prevents a small-cell-line model from being promoted on unstable performance.",
  "- This is a sensitivity analysis of in-vitro viability only. It does not establish clinical efficacy, direct target engagement, or an endothelial-specific mechanism."
), file.path(out, "drug_ml_report.md"))
