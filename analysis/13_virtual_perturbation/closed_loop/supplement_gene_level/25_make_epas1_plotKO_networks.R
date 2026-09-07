options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(scTenifoldKnk)
  library(igraph)
})

root <- "D:/PDAC_P1"
out <- file.path(root, "analysis/13_virtual_perturbation/closed_loop/supplement_gene_level")
datasets <- c("GSE154778", "GSE155698", "GSE212966")
results <- file.path(root, "analysis/13_virtual_perturbation/results")

draw_networks <- function(file, width = 12, height = 4.5) {
  if (grepl("\\.png$", file, ignore.case = TRUE)) {
    grDevices::png(file, width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::pdf(file, width = width, height = height)
  }
  old <- par(no.readonly = TRUE)
  on.exit({par(old); grDevices::dev.off()}, add = TRUE)
  par(mfrow = c(1, 3), mar = c(0.2, 0.2, 1.3, 0.2), oma = c(0, 0, 2.4, 0))
  for (dataset in datasets) {
    f <- file.path(results, paste0(dataset, "_primary__endothelial__EPAS1__KO.rds"))
    x <- readRDS(f)
    plotKO(x, "EPAS1", q = 0.99, annotate = FALSE, fdrThreshold = 0.05)
    title(main = dataset, line = 0.2, cex.main = 1.05)
    rm(x)
    gc(verbose = FALSE)
  }
  mtext("scTenifoldKnk EPAS1 knockout networks (model-inferred)", outer = TRUE,
        line = 0.9, cex = 1.2, font = 2)
}

draw_networks(file.path(out, "Figure_S27_EPAS1_plotKO_networks.pdf"))
draw_networks(file.path(out, "Figure_S27_EPAS1_plotKO_networks.png"), width = 12, height = 4.5)
