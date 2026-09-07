options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(edgeR); library(glmnet); library(pROC)})
root <- "D:/PDAC_P1"
out <- file.path(root, "analysis/14_prediction/endothelial_state")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
targets <- c("EPAS1", "TACC1", "MARCKS", "HERPUD1")
markers <- c("PECAM1", "VWF", "EMCN")

zscore_rows <- function(m) { m <- as.matrix(m); t(scale(t(m))) }
gene_matrix <- function(x, genes) {
  genes <- intersect(genes, rownames(x)); if (!length(genes)) return(matrix(numeric(),0,ncol(x)))
  x <- x[genes,,drop=FALSE]; x <- x[!duplicated(rownames(x)),,drop=FALSE]; x
}
make_state <- function(expr, dataset) {
  marker_mat <- gene_matrix(expr, markers)
  target_mat <- gene_matrix(expr, targets)
  if (nrow(marker_mat) < 2 || nrow(target_mat) < 2) stop(dataset, ": insufficient genes")
  marker_score <- colMeans(zscore_rows(marker_mat), na.rm=TRUE)
  q <- quantile(marker_score, c(.30,.70), na.rm=TRUE)
  label <- ifelse(marker_score >= q[2], 1L, ifelse(marker_score <= q[1], 0L, NA_integer_))
  list(dataset=dataset, expr=target_mat, marker_score=marker_score, label=label, thresholds=q)
}

# TCGA primary tumors: counts -> logCPM, map Ensembl annotation to symbols.
counts <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_primary_tumor_counts_unstranded.rds"))
ann <- readRDS(file.path(root,"analysis/01_tcga_qc/TCGA-PAAD_gene_annotation.rds"))
sym <- ann$gene_name[match(rownames(counts), ann$gene_id)]
keep <- !is.na(sym) & sym != ""; counts <- counts[keep,,drop=FALSE]; rownames(counts) <- make.unique(sym[keep])
tcga <- cpm(DGEList(counts=counts), log=TRUE, prior.count=2)
sets <- list(TCGA_PAAD=make_state(tcga,"TCGA-PAAD"))

# GSE71729 GEO-provided normalized matrix; rows are re-annotated symbols.
geo <- file.path(root,"data/01_bulk/GSE71729/raw/GSE71729_series_matrix.txt.gz")
ln <- readLines(gzfile(geo)); begin <- grep("^!series_matrix_table_begin",ln)[1]
dat <- fread(geo, skip=begin-1, data.table=TRUE); setnames(dat,1,"gene")
gexpr <- as.matrix(dat[,-1,with=FALSE]); rownames(gexpr) <- as.character(dat$gene)
gexpr <- gexpr[!duplicated(rownames(gexpr)),,drop=FALSE]
sets$GSE71729 <- make_state(gexpr,"GSE71729")

# GSE62452 RMA gene-symbol matrix.
ge624 <- fread(file.path(root,"analysis/08_program_validation/GSE62452/GSE62452_RMA_gene_symbol.tsv"), data.table=FALSE, check.names=FALSE)
rownames(ge624) <- ge624[[1]]; ge624 <- as.matrix(ge624[,-1,drop=FALSE])
sets$GSE62452 <- make_state(ge624,"GSE62452")

train <- sets$TCGA_PAAD; ix <- which(!is.na(train$label)); X <- t(train$expr[,ix,drop=FALSE]); y <- train$label[ix]
set.seed(20260903)
fit <- cv.glmnet(X,y,family="binomial",alpha=.5,nfolds=10,type.measure="auc",standardize=TRUE)
oof <- as.numeric(predict(fit,newx=X,s="lambda.1se",type="response"))
roc_train <- roc(y,oof,quiet=TRUE)
coef_tab <- as.matrix(coef(fit,s="lambda.1se")); coef_df <- data.frame(feature=rownames(coef_tab), coefficient=as.numeric(coef_tab[,1]))
fwrite(coef_df,file.path(out,"elastic_net_coefficients.tsv"),sep="\t")
pred_rows <- list(data.frame(dataset="TCGA-PAAD", sample=rownames(X), label=y, marker_score=train$marker_score[ix], predicted_probability=oof))
metric <- function(ds, truth, prob) { r <- roc(truth,prob,quiet=TRUE); data.frame(dataset=ds,n=length(truth),n_high=sum(truth==1),n_low=sum(truth==0),AUROC=as.numeric(auc(r)),AUPRC=as.numeric(pr.curve <- NA), stringsAsFactors=FALSE) }
metric_one <- function(ds, truth, prob) {
  r <- roc(truth, prob, quiet=TRUE)
  cal <- glm(truth ~ qlogis(pmin(pmax(prob,1e-6),1-1e-6)), family=binomial())
  data.frame(dataset=ds, n=length(truth), n_high=sum(truth==1), n_low=sum(truth==0),
             AUROC=as.numeric(auc(r)), AUROC_CI_low=ci.auc(r,boot.n=500,quiet=TRUE)[1],
             AUROC_CI_high=ci.auc(r,boot.n=500,quiet=TRUE)[3],
             Brier=mean((prob-truth)^2), calibration_intercept=unname(coef(glm(truth~1+offset(qlogis(pmin(pmax(prob,1e-6),1-1e-6))),family=binomial()))),
             calibration_slope=unname(coef(cal)[2]), stringsAsFactors=FALSE)
}
metrics <- list(metric_one("TCGA-PAAD", y, oof))
for (nm in setdiff(names(sets),"TCGA_PAAD")) {
  s <- sets[[nm]]; valid <- which(!is.na(s$label)); xx <- t(s$expr[,valid,drop=FALSE]); pp <- as.numeric(predict(fit,newx=xx,s="lambda.1se",type="response")); rr <- roc(s$label[valid],pp,quiet=TRUE)
  pred_rows[[length(pred_rows)+1]] <- data.frame(dataset=nm,sample=rownames(xx),label=s$label[valid],marker_score=s$marker_score[valid],predicted_probability=pp)
  metrics[[length(metrics)+1]] <- metric_one(nm, s$label[valid], pp)
}
pred <- rbindlist(pred_rows,fill=TRUE); fwrite(pred,file.path(out,"predictions.tsv"),sep="\t")
met <- rbindlist(metrics,fill=TRUE); fwrite(met,file.path(out,"performance.tsv"),sep="\t")

# Exact SHAP values for the final standardized elastic-net linear model.
# For a linear log-odds model, phi_j = beta_j*(x_j - mean_j); this is mathematically exact.
mu <- colMeans(X); sx <- apply(X,2,sd); beta <- coef(fit,s="lambda.1se")[-1,1]; names(beta) <- colnames(X)
Xstd <- sweep(sweep(X,2,mu,"-"),2,sx,"/"); shap <- sweep(Xstd,2,beta,"*")
shap_df <- data.frame(dataset="TCGA-PAAD",sample=rownames(X),shap,check.names=FALSE)
for (nm in setdiff(names(sets),"TCGA_PAAD")) {s<-sets[[nm]]; valid<-which(!is.na(s$label)); xx<-t(s$expr[,valid,drop=FALSE]); xx<-xx[,colnames(X),drop=FALSE]; xstd<-sweep(sweep(xx,2,mu,"-"),2,sx,"/"); sh<-sweep(xstd,2,beta,"*"); shap_df<-rbind(shap_df,data.frame(dataset=nm,sample=rownames(xx),sh,check.names=FALSE))}
fwrite(shap_df,file.path(out,"exact_linear_shap_values.tsv"),sep="\t")
writeLines(c(
  "# Endothelial-state prediction and SHAP audit", "",
  "## Design", "- Label: top/bottom 30% of an independent endothelial marker score (PECAM1/VWF/EMCN) within each cohort; middle 40% excluded.",
  "- Predictors: pre-specified converged candidates EPAS1, TACC1, MARCKS and HERPUD1.",
  "- Model: 10-fold cross-validated elastic-net logistic regression; TCGA-PAAD training with external GSE71729 and GSE62452 validation.",
  "- Interpretation: exact additive SHAP values for the final linear log-odds model; these are predictive contributions, not causal effects.", "",
  "## Gate", "SHAP enters the main manuscript only if external AUROC is reproducibly >=0.70 in both cohorts; otherwise retain the model and SHAP as exploratory supplement or omit it.", "",
  "## Outputs", "- performance.tsv", "- elastic_net_coefficients.tsv", "- predictions.tsv", "- exact_linear_shap_values.tsv"
), file.path(out,"prediction_shap_report.md"))
