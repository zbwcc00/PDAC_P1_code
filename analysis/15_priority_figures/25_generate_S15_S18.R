base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(patchwork)})
theme_pub <- theme_classic(base_size = 9, base_family = "Arial") + theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8.5, colour = "#555555"), axis.title = element_text(size = 9), axis.text = element_text(colour = "#222222"), legend.position = "bottom", legend.title = element_blank(), strip.background = element_rect(fill = "#F2F4F6", colour = NA), strip.text = element_text(face = "bold"))
save_pub <- function(p, stem, width, height) {ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height, device = cairo_pdf); png(file.path(out_dir, paste0(stem, ".png")), width = round(width * 300), height = round(height * 300), res = 300, type = "cairo-png"); print(p); dev.off()}
read_tsv <- function(x) read.delim(x, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")

dataset_files <- c(GSE154778 = "GSE154778__primary_tumor.rds", GSE155698 = "GSE155698__primary_tumor.rds", GSE212966 = "GSE212966__primary_tumor.rds")
set.seed(20260905)
umap_list <- list(); expr_list <- list()
for (ds in names(dataset_files)) {
  obj <- readRDS(file.path(base_dir, "analysis/03_scrna_annotation/seurat_annotated", dataset_files[[ds]]))
  obj <- subset(obj, cells = sample(colnames(obj), min(2500, ncol(obj))))
  if (length(VariableFeatures(obj)) == 0) obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 20, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
  emb <- as.data.frame(Embeddings(obj, "umap")); emb$major_lineage <- obj$major_lineage; emb$dataset <- ds; umap_list[[ds]] <- emb
  ex <- FetchData(obj, vars = c("EPAS1", "major_lineage")); ex$dataset <- ds; expr_list[[ds]] <- ex
}
umap_df <- do.call(rbind, umap_list); expr_df <- do.call(rbind, expr_list)
umap_df <- umap_df[!is.na(umap_df$major_lineage), ]
lineage_cols <- c(endothelial = "#0072B2", fibroblast_CAF = "#56B4E9", myeloid = "#009E73", T_NK = "#CC79A7", B_cell = "#E69F00", epithelial = "#999999", malignant_epithelial = "#264653", acinar = "#F0E442", mast = "#D55E00", unknown_ambiguous = "#BBBBBB")
umap_df$major_lineage <- factor(umap_df$major_lineage, levels = names(lineage_cols))
umap_df <- umap_df[!is.na(umap_df$major_lineage), ]
p15a <- ggplot(umap_df, aes(umap_1, umap_2, colour = major_lineage)) + geom_point(size = 0.22, alpha = 0.55) + facet_wrap(~dataset, nrow = 1) + scale_colour_manual(values = lineage_cols, drop = FALSE) + labs(title = "Annotated primary-tumor single-cell landscapes", x = "UMAP 1", y = "UMAP 2") + theme_pub + theme(legend.position = "right", legend.text = element_text(size = 7))
expr_df$major_lineage <- factor(expr_df$major_lineage, levels = names(lineage_cols))
expr_df <- expr_df[expr_df$major_lineage %in% c("endothelial", "malignant_epithelial", "fibroblast_CAF", "myeloid", "T_NK"), ]
expr_df$major_lineage <- factor(expr_df$major_lineage, levels = c("endothelial", "malignant_epithelial", "fibroblast_CAF", "myeloid", "T_NK"))
p15b <- ggplot(expr_df, aes(major_lineage, EPAS1, fill = major_lineage)) + geom_violin(scale = "width", trim = TRUE, colour = NA, alpha = 0.72) + geom_boxplot(width = 0.12, outlier.shape = NA, colour = "#333333") + facet_wrap(~dataset, nrow = 1) + scale_fill_manual(values = lineage_cols, drop = FALSE) + labs(title = "EPAS1 expression in major lineages", x = NULL, y = "EPAS1 normalized expression") + theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5), legend.position = "none")
cnv <- read_tsv(file.path(base_dir, "analysis/04_cnv_copykat/PDAC_copykat_all_sample_summary.tsv")); cnv$dataset <- gsub('"', "", cnv$dataset); cnv$patient_id <- gsub('"', "", cnv$patient_id); cnv$dataset <- factor(cnv$dataset, levels = names(dataset_files))
p15c <- ggplot(cnv, aes(dataset, aneuploid_fraction, colour = dataset)) + geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.65) + geom_boxplot(outlier.shape = NA, width = 0.42, colour = "#333333", fill = NA) + scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.08))) + scale_colour_manual(values = c(GSE154778 = "#0072B2", GSE155698 = "#009E73", GSE212966 = "#D55E00")) + labs(title = "CopyKAT aneuploid fraction", x = NULL, y = "Aneuploid cells / sampled cells") + theme_pub + theme(legend.position = "none")
run <- read_tsv(file.path(base_dir, "analysis/04_cnv_copykat/PDAC_copykat_run_summary.tsv")); run$dataset <- gsub('"', "", run$dataset)
p15d <- ggplot(run, aes(x = reorder(dataset, aneuploid_fraction), y = aneuploid_fraction, fill = dataset)) + geom_col(width = 0.62) + geom_text(aes(label = sprintf("%.1f%%", 100 * aneuploid_fraction)), vjust = -0.35, size = 3) + scale_y_continuous(limits = c(0, 0.52), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.10))) + scale_fill_manual(values = c(GSE154778 = "#0072B2", GSE155698 = "#009E73", GSE212966 = "#D55E00")) + labs(title = "Cohort-level CopyKAT summary", x = NULL, y = "Aneuploid fraction") + theme_pub + theme(legend.position = "none")
p15 <- (p15a | p15b) / (p15c | p15d) + plot_annotation(title = "Supplementary Figure S15. Single-cell annotation and CopyKAT/CNV assessment", subtitle = "UMAPs use 2,500 randomly sampled primary-tumor cells per cohort; CopyKAT panels retain patient/sample as the summary unit.") & theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.5))
save_pub(p15, "Figure_S15_scRNA_annotation_CopyKAT_CNV", 9.2, 7.2)

dock <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_exhaustiveness16_summary.tsv"))
dock$compound_label <- ifelse(dock$compound == "BRD-K56751279", "Y-39983", "triclabendazole")
dock$compound_label <- factor(dock$compound_label, levels = c("Y-39983", "triclabendazole"))
p16a <- ggplot(dock, aes(pdb_id, affinity_best_kcal_mol, colour = compound_label, group = compound_label)) + geom_line(linewidth = 0.7) + geom_point(size = 2.5) + scale_colour_manual(values = c("Y-39983" = "#D55E00", "triclabendazole" = "#0072B2")) + labs(title = "Five-conformer docking at exhaustiveness 16", x = "EPAS1 conformation", y = "Best-mode affinity (kcal/mol)") + theme_pub
contact <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_contact_by_conformer.tsv"))
contact$group <- factor(contact$group, levels = c("Y-39983", "triclabendazole", "negative_control"))
contact <- contact[!is.na(contact$group), ]
p16b <- ggplot(contact, aes(pdb_id, contact_count, fill = group)) + geom_col(position = position_dodge(width = 0.75), width = 0.65) + scale_fill_manual(values = c("Y-39983" = "#D55E00", "triclabendazole" = "#0072B2", "negative_control" = "#999999"), labels = c("Y-39983", "triclabendazole", "Negative control")) + labs(title = "Pose contact counts", x = "EPAS1 conformation", y = "Protein-ligand contacts") + theme_pub
cons <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_consensus_binding_residues.tsv"))
cons$group <- factor(cons$group, levels = c("Y-39983", "triclabendazole", "negative_control", "positive_control")); cons <- cons[!is.na(cons$group), ]; cons$residue <- sub("A:", "", cons$residue)
p16c <- ggplot(cons, aes(residue, group, fill = frequency)) + geom_tile(colour = "white") + scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", limits = c(0, 1)) + labs(title = "Consensus contact residue frequency", x = "EPAS1 residue", y = NULL, fill = "Across conformers") + theme_pub + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.text.y = element_text(size = 8))
qc <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/multiconformer_summary_corrected.tsv")); qc <- qc[qc$kind == "cocrystal_redock", ]
p16d <- ggplot(qc, aes(pdb_id, rmsd_angstrom)) + geom_hline(yintercept = 2, linetype = 2, colour = "#D55E00") + geom_col(fill = "#264653", width = 0.62) + geom_text(aes(label = sprintf("%.2f", rmsd_angstrom)), vjust = -0.35, size = 3) + scale_y_continuous(limits = c(0, 2.35), expand = expansion(mult = c(0, 0.08))) + labs(title = "Co-crystal redocking QC", x = "EPAS1 conformation", y = "Heavy-atom RMSD (A)") + theme_pub
p16 <- (p16a | p16b) / (p16c | p16d) + plot_annotation(title = "Supplementary Figure S16. EPAS1 structure evidence and docking QC", subtitle = "Candidate contact patterns are structural hypotheses; the dashed line marks the prespecified RMSD <=2 A pass threshold.") & theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.5))
save_pub(p16, "Figure_S16_EPAS1_structure_docking_QC", 8.6, 7.3)

mr <- read_tsv(file.path(base_dir, "analysis/11_methodology_review/supplemental_gate/MR_instrument_gate_summary.tsv"))
mr$gene <- factor(mr$gene, levels = rev(mr$gene)); mr$GTEx_overlap_n <- as.numeric(mr$GTEx_overlap_n); mr$BLUEPRINT_cis_n <- as.numeric(mr$BLUEPRINT_cis_n); mr$LEPIK_cis_n <- as.numeric(mr$LEPIK_cis_n)
mr_long <- rbind(data.frame(gene = mr$gene, resource = "GTEx v8", n = mr$GTEx_overlap_n), data.frame(gene = mr$gene, resource = "BLUEPRINT", n = mr$BLUEPRINT_cis_n), data.frame(gene = mr$gene, resource = "Lepik", n = mr$LEPIK_cis_n))
p17a <- ggplot(mr_long, aes(resource, gene, fill = n)) + geom_tile(colour = "white") + geom_text(aes(label = n), size = 3) + scale_fill_gradient(low = "#F7FBFF", high = "#2166AC") + labs(title = "cis-eQTL instrument availability", x = NULL, y = NULL, fill = "cis-eQTL count") + theme_pub
mr_status <- data.frame(gate = c("GTEx genome-wide instrument", "GTEx outcome overlap", "BLUEPRINT/Lepik", "Primary MR promotion"), status = c("No", "No valid set", "Exploratory only", "Failed"), detail = c("EPAS1: false", "EPAS1: n=0", "Tissue mismatch / limited overlap", "No causal claim")); mr_status$gate <- factor(mr_status$gate, levels = rev(mr_status$gate))
p17b <- ggplot(mr_status, aes(x = 1, y = gate, colour = status)) + geom_point(size = 5) + geom_text(aes(label = detail), hjust = -0.1, size = 3.1, colour = "#222222") + scale_x_continuous(limits = c(0.8, 2.4), breaks = NULL) + scale_colour_manual(values = c("No" = "#D55E00", "No valid set" = "#D55E00", "Exploratory only" = "#E69F00", "Failed" = "#D55E00")) + labs(title = "Promotion gate", x = NULL, y = NULL) + theme_pub + theme(legend.position = "none")
p17 <- (p17a | p17b) + plot_annotation(title = "Supplementary Figure S17. MR instrument and promotion-gate audit") & theme(plot.title = element_text(face = "bold", size = 12))
save_pub(p17, "Figure_S17_MR_instrument_gate", 8.2, 4.6)

surv <- read_tsv(file.path(base_dir, "analysis/08_program_validation/GSE21501/GSE21501_EPAS1_external_survival.tsv")); surv$feature_label <- ifelse(surv$feature == "EPAS1_expression", "EPAS1 expression", "Endothelial program"); surv$se <- as.numeric(surv$SE); surv$HR <- as.numeric(surv$HR); surv$low <- exp(log(surv$HR) - 1.96 * surv$se); surv$high <- exp(log(surv$HR) + 1.96 * surv$se); surv$feature_label <- factor(surv$feature_label, levels = rev(c("EPAS1 expression", "Endothelial program"))); surv$label_x <- 1.42
p18 <- ggplot(surv, aes(HR, feature_label)) + geom_vline(xintercept = 1, linetype = 2, colour = "#777777") + geom_errorbarh(aes(xmin = low, xmax = high), height = 0.2, linewidth = 0.9, colour = "#264653") + geom_point(size = 3, colour = "#D55E00") + geom_text(aes(x = label_x, label = sprintf("HR=%.3f; P=%.3f", HR, P)), hjust = 0, size = 3) + scale_x_log10(limits = c(0.35, 2.3), breaks = c(0.5, 1, 2)) + labs(title = "Supplementary Figure S18. GSE21501 external survival sensitivity", x = "Hazard ratio (log scale)", y = NULL, subtitle = "n=102 samples; 66 events; neither prespecified association was significant") + theme_pub
save_pub(p18, "Figure_S18_GSE21501_survival_forest", 7.0, 3.8)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "25_R_sessionInfo.txt"))
