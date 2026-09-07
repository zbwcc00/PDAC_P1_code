options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

outdir <- 'D:/PDAC_P1/analysis/13_virtual_perturbation'
objdir <- 'D:/PDAC_P1/analysis/03_scrna_annotation/seurat_annotated'
files <- c(
  GSE154778_primary='GSE154778__primary_tumor.rds',
  GSE155698_primary='GSE155698__primary_tumor.rds',
  GSE212966_primary='GSE212966__primary_tumor.rds'
)
candidates <- c('EPAS1','TACC1','MARCKS','HERPUD1','CTHRC1')

find_celltype <- function(md) {
  candidates <- c('cell_type','celltype','CellType','celltype_major','major_celltype','major_lineage',
                  'annotation','annotated_cell_type','ident','seurat_annotations')
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit)) return(hit[1])
  return(NA_character_)
}

all_cov <- list(); all_type <- list(); i <- 0L
for (nm in names(files)) {
  f <- file.path(objdir, files[[nm]])
  message('Reading ', nm)
  obj <- readRDS(f)
  assay <- if ('RNA' %in% names(obj@assays)) 'RNA' else DefaultAssay(obj)
  genes <- rownames(obj[[assay]])
  md <- obj@meta.data
  ctcol <- find_celltype(md)
  if (is.na(ctcol)) {
    md$.__celltype__ <- as.character(Idents(obj)); ctcol <- '.__celltype__'
  }
  ct <- as.character(md[[ctcol]])
  names(ct) <- rownames(md)
  dat <- GetAssayData(obj, assay=assay, layer='data')
  present <- intersect(candidates, rownames(dat))
  for (g in candidates) {
    i <- i + 1L
    if (!(g %in% rownames(dat))) {
      all_cov[[i]] <- data.frame(dataset=nm, assay=assay, gene=g, n_cells=ncol(obj), n_detected=0,
                                  detection_rate=0, n_types=0, celltype_column=ctcol)
      next
    }
    x <- as.numeric(dat[g, ])
    all_cov[[i]] <- data.frame(dataset=nm, assay=assay, gene=g, n_cells=length(x),
                                n_detected=sum(x > 0), detection_rate=mean(x > 0),
                                n_types=length(unique(ct[which(x > 0)])), celltype_column=ctcol)
    tab <- aggregate(x > 0, by=list(cell_type=ct), FUN=sum)
    names(tab)[2] <- 'n_detected'
    denom <- as.data.frame(table(ct), stringsAsFactors=FALSE)
    names(denom) <- c('cell_type','n_cells')
    tab <- merge(denom, tab, by='cell_type', all.x=TRUE)
    tab$n_detected[is.na(tab$n_detected)] <- 0
    tab$detection_rate <- tab$n_detected / tab$n_cells
    tab$dataset <- nm; tab$gene <- g; tab$assay <- assay; tab$celltype_column <- ctcol
    all_type[[length(all_type)+1L]] <- tab[,c('dataset','assay','gene','cell_type','n_cells','n_detected','detection_rate','celltype_column')]
  }
  rm(obj, dat); gc()
}
cov <- do.call(rbind, all_cov)
typ <- if (length(all_type)) do.call(rbind, all_type) else data.frame()
write.table(cov, file.path(outdir,'candidate_gene_coverage.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(typ, file.path(outdir,'candidate_gene_coverage_by_celltype.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir,'coverage_sessionInfo.txt'))
print(cov)
