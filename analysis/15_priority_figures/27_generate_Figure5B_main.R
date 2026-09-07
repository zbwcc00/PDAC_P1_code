options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

base_dir <- "D:/PDAC_P1"
out_dir <- file.path(base_dir, "analysis/15_priority_figures")
read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")
theme_pub <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, colour = "#555555"),
        axis.title = element_text(size = 9), axis.text = element_text(colour = "#222222"),
        legend.position = "bottom", legend.title = element_blank())

dock <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_exhaustiveness16_summary.tsv"))
dock$compound_label <- factor(ifelse(dock$compound == "BRD-K56751279", "Y-39983", "triclabendazole"),
                              levels = c("Y-39983", "triclabendazole"))
p_affinity <- ggplot(dock, aes(pdb_id, affinity_best_kcal_mol, colour = compound_label, group = compound_label)) +
  geom_line(linewidth = 0.75) + geom_point(size = 2.5) +
  scale_colour_manual(values = c("Y-39983" = "#D55E00", "triclabendazole" = "#0072B2")) +
  labs(title = "Five-conformer Vina affinity", x = "EPAS1 conformation", y = "Best-mode affinity (kcal/mol)") + theme_pub

qc <- read_tsv(file.path(base_dir, "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/multiconformer_summary_corrected.tsv"))
qc <- qc[qc$kind == "cocrystal_redock", ]
qc$rmsd_angstrom <- as.numeric(qc$rmsd_angstrom)
p_rmsd <- ggplot(qc, aes(pdb_id, rmsd_angstrom)) +
  geom_hline(yintercept = 2, linetype = 2, colour = "#D55E00", linewidth = .55) +
  geom_col(fill = "#264653", width = .62) +
  geom_text(aes(label = sprintf("%.2f", rmsd_angstrom)), vjust = -0.35, size = 3) +
  scale_y_continuous(limits = c(0, 2.35), expand = expansion(mult = c(0, .08))) +
  labs(title = "Co-crystal redocking QC", x = "EPAS1 conformation", y = "Heavy-atom RMSD (A)") + theme_pub

p <- (p_affinity | p_rmsd) + plot_annotation(
  title = "Figure 5B. Multi-conformer docking stability and redocking quality control",
  subtitle = "Five EPAS1 PAS-B conformations; dashed line indicates the prespecified RMSD <= 2 A pass threshold.",
  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 14),
                plot.subtitle = element_text(family = "Arial", size = 9, colour = "#555555")))

stem <- file.path(out_dir, "Figure_S16_EPAS1_main_docking_QC")
ggsave(paste0(stem, ".pdf"), p, width = 8.6, height = 3.8, device = cairo_pdf)
png(paste0(stem, ".png"), width = 2580, height = 1140, res = 300, type = "cairo-png")
print(p)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out_dir, "27_Figure5B_R_sessionInfo.txt"))
