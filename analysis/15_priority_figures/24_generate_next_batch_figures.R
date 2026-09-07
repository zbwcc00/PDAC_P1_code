base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "#555555"),
        axis.title = element_text(size = 9),
        axis.text = element_text(colour = "#222222"),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "#F2F4F6", colour = NA),
        strip.text = element_text(face = "bold"))

save_pub <- function(plot_obj, stem, width = 7.2, height = 4.6) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot_obj, width = width, height = height, device = cairo_pdf)
  png(file.path(out_dir, paste0(stem, ".png")), width = round(width * 300), height = round(height * 300), res = 300, type = "cairo-png")
  print(plot_obj)
  dev.off()
}

read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")

protein_detection <- read_tsv(file.path(base_dir, "analysis/11_methodology_review/CPTAC_PAAD/CPTAC_protein_detection_summary.tsv"))
protein_cox <- read_tsv(file.path(base_dir, "analysis/11_methodology_review/CPTAC_PAAD/CPTAC_protein_survival_cox.tsv"))
protein_cox$gene <- gsub('"', "", protein_cox$gene)
protein_detection$gene <- factor(protein_detection$gene, levels = rev(protein_detection$gene))
protein_detection$IQR_low <- protein_detection$median - protein_detection$IQR / 2
protein_detection$IQR_high <- protein_detection$median + protein_detection$IQR / 2
protein_cox$gene <- factor(protein_cox$gene, levels = rev(protein_cox$gene))
protein_cox$se <- abs(log(protein_cox$continuous_HR)) / qnorm(1 - protein_cox$continuous_P / 2)
protein_cox$low <- exp(log(protein_cox$continuous_HR) - 1.96 * protein_cox$se)
protein_cox$high <- exp(log(protein_cox$continuous_HR) + 1.96 * protein_cox$se)

p11a <- ggplot(protein_detection, aes(x = median, y = gene)) +
  geom_errorbarh(aes(xmin = IQR_low, xmax = IQR_high), height = 0.22, colour = "#0072B2", linewidth = 0.8) +
  geom_point(size = 3, colour = "#D55E00") +
  geom_text(aes(label = sprintf("%.2f", median)), hjust = -0.25, size = 3.0) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(title = "CPTAC protein detection", x = "Median protein abundance (log2 scale)", y = NULL,
       subtitle = "n = 140 samples; bars denote median ± IQR/2") + theme_pub

p11b <- ggplot(protein_cox, aes(y = gene, x = continuous_HR)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.20, colour = "#264653", linewidth = 0.8) +
  geom_point(size = 3, colour = "#264653") +
  geom_point(aes(x = median_split_HR), shape = 2, size = 2.8, colour = "#D55E00") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2), limits = c(0.2, 2.2)) +
  labs(title = "Cox sensitivity analysis", x = "Hazard ratio (log scale)", y = NULL,
       subtitle = "Filled: continuous; open triangle: median split") + theme_pub

p11 <- p11a + p11b + plot_annotation(title = "Supplementary Figure S11. CPTAC protein-layer evidence",
                                      subtitle = "Protein availability was strong; prognostic inference remained sensitivity-dependent.") &
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9))
save_pub(p11, "Figure_S11_CPTAC_protein_evidence", 7.4, 4.5)

priority <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/drugreflector/drugreflector_priority.tsv"))
mapping <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/drugreflector/drugreflector_compound_mapping_pubchem.tsv"))
priority <- priority[priority$target == "EPAS1", ]
priority <- priority[order(priority$mean_rank), ]
priority <- priority[seq_len(min(15, nrow(priority))), ]
priority$name <- priority$compound
idx <- match(priority$compound, mapping$compound)
priority$name[!is.na(idx) & nzchar(mapping$cmap_name[idx])] <- mapping$cmap_name[idx][!is.na(idx) & nzchar(mapping$cmap_name[idx])]
priority$name[priority$compound == "BRD-K56751279"] <- "Y-39983"
priority$name[priority$compound == "BRD-K81916719"] <- "triclabendazole"
priority$highlight <- ifelse(priority$name %in% c("Y-39983", "triclabendazole"), "Prespecified docking candidate", "Other recurrent hit")
priority$name <- factor(priority$name, levels = rev(priority$name))

p12 <- ggplot(priority, aes(x = mean_probability, y = name, fill = highlight)) +
  geom_col(width = 0.68, colour = "white") +
  geom_text(aes(label = sprintf("rank %d", round(mean_rank) + 1)), hjust = -0.12, size = 3.0) +
  scale_fill_manual(values = c("Prespecified docking candidate" = "#D55E00", "Other recurrent hit" = "#8C8C8C")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "EPAS1 reverse-signature prioritization", x = "Mean reversal probability", y = NULL,
       subtitle = "Top 15 hits; each displayed hit recurred in both prespecified strata") + theme_pub
save_pub(p12, "Figure_S12_DrugReflector_prioritization", 7.2, 5.4)

dep <- read_tsv(file.path(base_dir, "analysis/11_methodology_review/DepMap_PDAC_candidate_dependency_summary.tsv"))
dep <- dep[dep$gene_id %in% c("EPAS1", "MARCKS", "HERPUD1", "TACC1"), ]
dep$gene_id <- factor(dep$gene_id, levels = rev(dep$gene_id))
p13a <- ggplot(dep, aes(x = median_gene_effect, y = gene_id)) +
  geom_vline(xintercept = c(-1, -0.5), linetype = 2, colour = "#999999") +
  geom_point(aes(size = frac_dependency_leq_neg05), colour = "#0072B2") +
  scale_size_continuous(range = c(2.5, 6), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "DepMap PDAC dependency", x = "Median gene effect", y = NULL, size = "Fraction ≤ −0.5") + theme_pub

crispr <- read_tsv(file.path(base_dir, "analysis/11_methodology_review/DepMap_PRISM_pharmacology/DepMap_CRISPR_PRISM_association_results.tsv"))
prism <- read_tsv(file.path(base_dir, "analysis/11_methodology_review/DepMap_PRISM_pharmacology/DepMap_PRISM_association_results.tsv"))
prism <- prism[prism$feature %in% c("EPAS1", "MARCKS", "HERPUD1", "TACC1"), ]
prism$drug_label <- ifelse(grepl("triclabendazole", prism$drug), "triclabendazole", "Y-39983")
prism$gene <- factor(prism$feature, levels = rev(c("EPAS1", "MARCKS", "HERPUD1", "TACC1")))
p13b <- ggplot(prism, aes(x = drug_label, y = gene, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("rho=%.2f\nFDR=%.2f", spearman_rho, FDR)), size = 2.8, lineheight = 0.9) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-0.5, 0.5)) +
  labs(title = "PRISM association gate", x = NULL, y = NULL, fill = "Spearman ρ") + theme_pub

ml <- read_tsv(file.path(base_dir, "analysis/14_prediction/DepMap_PRISM_drug_ML/drug_ml_repeated_CV_summary.tsv"))
ml$drug_label <- ifelse(grepl("triclabendazole", ml$drug), "triclabendazole", "Y-39983")
p13c <- ggplot(ml, aes(x = median_elastic_net_R2, y = drug_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = R2_2.5pct, xmax = R2_97.5pct), height = 0.20, colour = "#D55E00", linewidth = 0.9) +
  geom_point(size = 3, colour = "#D55E00") +
  labs(title = "Repeated-CV drug ML", x = "Elastic-net R² (median; 95% repetition interval)", y = NULL) + theme_pub

stability <- read_tsv(file.path(base_dir, "analysis/14_prediction/DepMap_PRISM_drug_ML/drug_ml_feature_stability.tsv"))
stability$drug_label <- ifelse(grepl("triclabendazole", stability$drug), "triclabendazole", "Y-39983")
stability <- stability[stability$feature %in% c("EPAS1", "MARCKS", "HERPUD1", "TACC1"), ]
stability$feature <- factor(stability$feature, levels = rev(c("EPAS1", "MARCKS", "HERPUD1", "TACC1")))
p13d <- ggplot(stability, aes(x = selection_frequency, y = feature, fill = drug_label)) +
  geom_point(aes(colour = drug_label), position = position_dodge(width = 0.7), size = 2.8) +
  scale_fill_manual(values = c("triclabendazole" = "#56B4E9", "Y-39983" = "#009E73")) +
  scale_colour_manual(values = c("triclabendazole" = "#56B4E9", "Y-39983" = "#009E73")) +
  scale_x_continuous(limits = c(-0.02, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Feature-selection stability", x = "Selection frequency", y = NULL) + theme_pub

p13 <- (p13a | p13b) / (p13c | p13d) + plot_annotation(title = "Supplementary Figure S13. Pharmacology and prediction negative gates",
                                                          subtitle = "These analyses define evidence boundaries and do not support efficacy or dependency claims.") &
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9))
save_pub(p13, "Figure_S13_negative_pharmacology_gates", 8.2, 7.1)

surv <- read_tsv(file.path(base_dir, "analysis/14_prediction/epas1_survival_ml/survival_model_performance.tsv"))
surv_long <- rbind(data.frame(model = surv$model, split = "TCGA internal", cindex = surv$TCGA_5fold_Cindex),
                   data.frame(model = surv$model, split = "CPTAC external", cindex = surv$CPTAC_external_Cindex))
surv_long$model <- factor(surv_long$model, levels = c("clinical", "molecular_cox", "elastic_net", "survival_forest"))
levels(surv_long$model) <- c("Clinical", "Molecular Cox", "Elastic-net", "Survival forest")
p14 <- ggplot(surv_long, aes(x = model, y = cindex, fill = split)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "#777777") +
  geom_col(position = position_dodge(width = 0.75), width = 0.62, colour = "white") +
  geom_text(aes(label = sprintf("%.3f", cindex)), position = position_dodge(width = 0.75), vjust = -0.35, size = 3) +
  scale_fill_manual(values = c("TCGA internal" = "#0072B2", "CPTAC external" = "#D55E00")) +
  scale_y_continuous(breaks = seq(0.4, 0.6, 0.05)) + coord_cartesian(ylim = c(0.4, 0.58)) +
  labs(title = "Survival-ML sensitivity analysis", x = NULL, y = "Harrell's C-index", fill = NULL,
       subtitle = "External discrimination remained near chance") + theme_pub +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_pub(p14, "Figure_S14_survival_ML_sensitivity", 7.0, 4.6)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "24_R_sessionInfo.txt"))

scrna <- read_tsv(file.path(base_dir, "analysis/02_scrna_qc/PDAC_scRNA_dataset_summary.tsv"))
comp <- read_tsv(file.path(base_dir, "analysis/03_scrna_annotation/PDAC_scRNA_patient_lineage_summary.tsv"))
scrna$dataset <- gsub('"', "", scrna$dataset)
comp$dataset <- gsub('"', "", comp$dataset)
scrna$dataset <- factor(scrna$dataset, levels = scrna$dataset)
lineages <- c("endothelial", "fibroblast_CAF", "myeloid", "T_NK", "B_cell", "epithelial", "malignant_epithelial", "acinar", "mast", "unknown_ambiguous")
mean_comp <- do.call(rbind, lapply(levels(scrna$dataset), function(ds) colMeans(comp[comp$dataset == ds, paste0("prop_", lineages)], na.rm = TRUE)))
mean_comp <- as.data.frame(mean_comp)
mean_comp$dataset <- factor(levels(scrna$dataset), levels = levels(scrna$dataset))
scale_factor <- max(scrna$samples) / max(scrna$total_cells / 1000)
p6a <- ggplot(scrna, aes(x = dataset, y = samples)) +
  geom_col(fill = "#0072B2", width = 0.68) +
  geom_line(aes(y = total_cells / 1000 * scale_factor, group = 1), colour = "#D55E00", linewidth = 1.1) +
  geom_point(aes(y = total_cells / 1000 * scale_factor), colour = "#D55E00", size = 2.8) +
  geom_text(aes(label = paste0("n=", samples)), vjust = -0.35, size = 3.0) +
  scale_y_continuous(name = "Samples", expand = expansion(mult = c(0, 0.12)), sec.axis = sec_axis(~ . / scale_factor, name = "Cells (x1000)")) +
  labs(title = "Discovery cohort scale", x = NULL) + theme_pub + theme(axis.text.x = element_text(angle = 25, hjust = 1))
mean_long <- reshape(mean_comp, varying = paste0("prop_", lineages), v.names = "fraction", timevar = "lineage", times = lineages, direction = "long")
mean_long$lineage <- factor(mean_long$lineage, levels = lineages)
p6b <- ggplot(mean_long, aes(x = dataset, y = fraction, fill = lineage)) +
  geom_col(width = 0.68) +
  scale_fill_manual(values = c("#0072B2", "#56B4E9", "#009E73", "#CC79A7", "#E69F00", "#999999", "#264653", "#F0E442", "#D55E00", "#BBBBBB"), labels = lineages) +
  labs(title = "Cellular composition", x = NULL, y = "Mean patient cell fraction") + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + theme_pub + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right", legend.text = element_text(size = 7))
p6 <- p6a + p6b + plot_annotation(title = "Supplementary Figure S6. Multi-cohort single-cell discovery overview") & theme(plot.title = element_text(face = "bold", size = 12))
save_pub(p6, "Figure_S6_scRNA_discovery_overview", 7.4, 4.2)
