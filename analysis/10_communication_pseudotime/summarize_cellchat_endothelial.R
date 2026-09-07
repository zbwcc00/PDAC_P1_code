options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "D:/PDAC_P1/analysis/10_communication_pseudotime/CellChat"
files_in <- list.files(root, pattern = "_to_endothelial\\.tsv$", full.names = TRUE)
files_out <- list.files(root, pattern = "_from_endothelial\\.tsv$", full.names = TRUE)
read_one <- function(f, direction) {x <- fread(f); x[, dataset := sub("_to_endothelial\\.tsv$|_from_endothelial\\.tsv$", "", basename(f))]; x[, direction := direction]; x}
all <- rbindlist(c(lapply(files_in, read_one, direction = "to_endothelial"), lapply(files_out, read_one, direction = "from_endothelial")), fill = TRUE)
all_sig <- all[pval < 0.05]
cons <- all_sig[, .(n_datasets = uniqueN(dataset), datasets = paste(sort(unique(dataset)), collapse = ","), mean_prob = mean(prob, na.rm=TRUE), min_p = min(pval, na.rm=TRUE)), by = .(direction, source, target, ligand, receptor, pathway_name)]
setorder(cons, -n_datasets, -mean_prob)
fwrite(cons, file.path(dirname(root), "CellChat_endothelial_consensus.tsv"), sep="\t")
top <- cons[n_datasets >= 2]
sink(file.path(dirname(root), "CellChat_endothelial_consensus_report.md")); cat("# CellChat endothelial consensus\\n\\n"); cat("Significant interactions (p<0.05) were retained per dataset; consensus requires recurrence in at least two datasets. This is exploratory and does not replace patient-level permutation tests.\\n\\n"); print(top); sink()
