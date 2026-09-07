options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
root <- "T:/"
out_dir <- file.path(root, "analysis/02_scrna_qc")
qc <- fread(file.path(out_dir, "PDAC_scRNA_cell_QC_summary.tsv"))
ctx <- fread(file.path(out_dir, "PDAC_scRNA_context_QC_summary.tsv"))
qc[, pass_fraction := cells_pass_core / cells_input]
qc[, review_flag := fifelse(median_genes < 300 | median_percent_mt > 30 | pass_fraction < 0.70, "exclude_or_manual_review",
                            fifelse(median_genes < 600 | pass_fraction < 0.85, "caution", "acceptable"))]
setorder(qc, review_flag, pass_fraction)
fwrite(qc, file.path(out_dir, "PDAC_scRNA_cell_QC_summary.tsv"), sep = "\t", na = "NA")
flags <- qc[review_flag != "acceptable", .(dataset, gsm, title, tissue_context, cells_input, cells_pass_core, pass_fraction, median_genes, median_percent_mt, doublet_rate, review_flag)]
fwrite(flags, file.path(out_dir, "PDAC_scRNA_sample_review_flags.tsv"), sep = "\t", na = "NA")

fmt <- function(x) format(round(x, 3), nsmall = 3, trim = TRUE)
lines <- c(
  "# PDAC P1 单细胞数据：元数据与 QC 门控报告",
  "",
  paste0("生成日期：", Sys.Date()),
  "",
  "## 1. 数据可用性",
  "",
  "- 统一 manifest 来自 GEO family SOFT 文件，保留全部公开样本，并将矩阵缺失与组织上下文显式标记。",
  "- GSE154778：16/16 个矩阵可读；10 个 primary tumor、6 个 metastasis。",
  "- GSE155698：38/41 个矩阵可读；15 个 primary tumor、3 个 adjacent normal、16 个 PDAC-PBMC、4 个 healthy-PBMC；缺失 GSM4710703、GSM4710705、GSM4710725 的可读矩阵。",
  "- GSE212966：12/12 个矩阵可读；6 个 PDAC、6 个 adjacent normal。",
  "",
  "## 2. QC 策略",
  "",
  "- 每个样本独立计算 detected genes、UMI 和 percent.mt；阈值由样本内 median/MAD 与 99th percentile 推导，而非跨数据集固定。",
  "- 核心过滤只使用 detected genes、UMI、percent.mt；高 UMI/高基因细胞仅作 review 信号。",
  "- `scDblFinder` 1.20.2 按样本运行；保留 doublet score/class，并将 `passes_core_QC` 写入 Seurat 元数据，不物理删除原始细胞。",
  "- 所有 Seurat 检查点保留 raw counts；对象按 dataset + tissue_context 分开保存，避免把 PBMC、adjacent 和 tumor 混在同一统计单元。",
  "",
  "## 3. 上下文汇总",
  "",
  "| Dataset | Context | Samples | Input cells | Core-pass cells | Median pass fraction | Median genes | Median doublet rate |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  paste0("| ", ctx$dataset, " | ", ctx$tissue_context, " | ", ctx$samples, " | ", ctx$cells_input, " | ", ctx$cells_pass_core, " | ", fmt(ctx$median_pass_fraction), " | ", fmt(ctx$median_genes), " | ", fmt(ctx$median_doublet_rate), " |"),
  "",
  "## 4. 必须人工复核的样本",
  "",
  if (nrow(flags)) paste0("- `", flags$dataset, "/", flags$gsm, "`（", flags$tissue_context, "）：median genes=", fmt(flags$median_genes), ", core-pass fraction=", fmt(flags$pass_fraction), ", median mt%=", fmt(flags$median_percent_mt), "；标记为 ", flags$review_flag, "。") else "- 未发现需要人工复核的样本。",
  "",
  "最重要的 QC 异常是 GSE154778 的 GSM4679533（P02）和 GSM4679544（MET03）：复杂度极低，应从发现性细胞组成与伪 bulk 分析中排除或单独做敏感性分析；不能仅凭细胞数量保留。",
  "",
  "## 5. P1 分析门控结论",
  "",
  "1. 发现集：GSE154778 primary tumor（排除 P02 后做主分析，含敏感性分析）与 GSE155698 primary tumor。",
  "2. 独立支持集：GSE212966 PDAC；其 adjacent 组仅作组织背景参考，不与 tumor 做细胞级统计混合。",
  "3. 转移生态支持：GSE154778 metastasis；用于检验恶性上皮程序与 CAF/髓系/CD8 生态是否在转移灶保持。",
  "4. PBMC 组：GSE155698 PDAC/healthy PBMC 只用于外周免疫背景，不作为肿瘤微环境 discovery 组。",
  "5. 下一步必须先完成 major lineage annotation、patient-level pseudobulk/比例汇总和预定义 CAF-myeloid-CD8 分数；若跨数据集方向一致，再进入 bulk 外部验证、空间定位、MR/共定位、虚拟扰动与药物反转。",
  "",
  "## 6. 解释边界",
  "",
  "- 该 QC 结果证明数据可进入细胞注释与生态关联分析，不等于证明免疫排斥机制、MR 因果关系或药物疗效。",
  "- 细胞数不是样本量；后续推断以患者/样本为统计单位，并严格区分 discovery、replication 和 treatment-context。"
)
writeLines(lines, file.path(out_dir, "PDAC_P1_scrna_metadata_qc_gate_report.md"), useBytes = TRUE)
