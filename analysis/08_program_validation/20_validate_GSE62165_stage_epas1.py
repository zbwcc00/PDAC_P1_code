from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr

BASE = Path(r"D:\PDAC_P1")
RAW = BASE / "data" / "01_bulk" / "GSE62165" / "raw"
OUT = BASE / "analysis" / "08_program_validation" / "GSE62165"
OUT.mkdir(parents=True, exist_ok=True)


def read_metadata(path: Path) -> pd.DataFrame:
    lines = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                lines.append(line.rstrip("\n"))
            elif line.startswith("!series_matrix_table_begin"):
                break
    fields = {}
    for line in lines:
        key, rest = line.split("\t", 1)
        values = rest.split("\t")
        fields.setdefault(key.lstrip("!"), []).append([v.strip('"') for v in values])
    accessions = fields["Sample_geo_accession"][0]
    meta = pd.DataFrame({"accession": accessions})
    chars = fields.get("Sample_characteristics_ch1", [])
    # GEO stores one row per characteristic field; recover stage/tissue by position.
    for row in chars:
        if len(row) != len(accessions):
            continue
        if any("grouped stage:" in x for x in row):
            meta["grouped_stage"] = [re.sub(r"^grouped stage:\s*", "", x, flags=re.I) for x in row]
        if any("tissue:" in x for x in row):
            meta["tissue"] = [re.sub(r"^tissue:\s*", "", x, flags=re.I) for x in row]
    title = fields.get("Sample_title", [[""] * len(meta)])[0]
    meta["title"] = title
    return meta


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip", comment="!", skiprows=68)


def main():
    matrix = read_matrix(RAW / "GSE62165_series_matrix_full.txt.gz")
    matrix = matrix.rename(columns={"ID_REF": "probe"})
    meta = read_metadata(RAW / "GSE62165_series_matrix_full.txt.gz")
    sample_cols = meta["accession"].tolist()
    matrix = matrix[["probe"] + sample_cols]
    ann_path = RAW / "GPL13667_full.txt"
    with ann_path.open("r", encoding="utf-8", errors="replace") as handle:
        ann_skip = next(i for i, line in enumerate(handle) if line.startswith("!platform_table_begin")) + 1
    ann = pd.read_csv(ann_path, sep="\t", skiprows=ann_skip, dtype=str)
    ann = ann.rename(columns={ann.columns[0]: "probe"})
    symbol_col = next((c for c in ann.columns if c.lower() == "gene symbol"), None)
    if symbol_col is None:
        symbol_col = next(c for c in ann.columns if "symbol" in c.lower())
    ann = ann[["probe", symbol_col]].rename(columns={symbol_col: "gene"})
    ann["gene"] = ann["gene"].fillna("").str.strip()
    matrix = matrix.merge(ann, on="probe", how="left")
    matrix = matrix[matrix["gene"].ne("") & matrix["gene"].ne("---")]
    gene_expr = matrix.groupby("gene", sort=False)[sample_cols].mean()

    cand = pd.read_csv(BASE / "analysis" / "07_candidate_screen" / "candidate_genes_stable_high_confidence_all.tsv", sep="\t")
    endo = cand[(cand["major_lineage"] == "endothelial") & (cand["stable_high_confidence"] == True)]
    up = [g for g in endo.loc[endo["direction"] == "tumor_up", "gene_id"] if g in gene_expr.index]
    down = [g for g in endo.loc[endo["direction"] == "tumor_down", "gene_id"] if g in gene_expr.index]
    def signed_score(genes_up, genes_down):
        z = gene_expr.loc[genes_up + genes_down].sub(gene_expr.loc[genes_up + genes_down].mean(axis=1), axis=0)
        z = z.div(gene_expr.loc[genes_up + genes_down].std(axis=1).replace(0, 1), axis=0)
        return z.loc[genes_up].mean(axis=0).sub(z.loc[genes_down].mean(axis=0), fill_value=0)
    meta["endothelial_score"] = signed_score(up, down).reindex(meta["accession"]).to_numpy()
    meta["EPAS1_expression"] = gene_expr.loc["EPAS1", meta["accession"]].to_numpy() if "EPAS1" in gene_expr.index else np.nan
    meta["is_tumor"] = meta["tissue"].str.contains("tumor", case=False, na=False) & ~meta["tissue"].str.contains("non-tumoral", case=False, na=False)
    meta.to_csv(OUT / "GSE62165_EPAS1_endothelial_stage_scores.tsv", sep="\t", index=False)

    rows = []
    valid_stage = meta[meta["grouped_stage"].notna() & ~meta["grouped_stage"].isin(["NA", "", "nan"])].copy()
    valid_stage = valid_stage[valid_stage["grouped_stage"].isin(["Early", "LNM", "Advanced"])]
    for variable in ["EPAS1_expression", "endothelial_score"]:
        groups = [g[variable].dropna().to_numpy() for _, g in valid_stage.groupby("grouped_stage")]
        labels = list(valid_stage.groupby("grouped_stage").groups)
        stat, p = kruskal(*groups) if len(groups) >= 2 else (np.nan, np.nan)
        rows.append({"analysis": "stage_kruskal", "variable": variable, "n": len(valid_stage), "groups": ";".join(labels), "statistic": stat, "p_value": p})
        for a, b in [("Early", "LNM"), ("Early", "Advanced"), ("LNM", "Advanced")]:
            x = valid_stage.loc[valid_stage["grouped_stage"] == a, variable].dropna()
            y = valid_stage.loc[valid_stage["grouped_stage"] == b, variable].dropna()
            if len(x) and len(y):
                u, p_pair = mannwhitneyu(x, y, alternative="two-sided")
                rows.append({"analysis": f"stage_pair_{a}_vs_{b}", "variable": variable, "n": len(x) + len(y), "groups": f"{a};{b}", "statistic": u, "p_value": p_pair})
    tumor = meta[meta["is_tumor"] & meta["EPAS1_expression"].notna()]
    if len(tumor) >= 3:
        rho, p = spearmanr(tumor["EPAS1_expression"], tumor["endothelial_score"])
        rows.append({"analysis": "tumor_gene_program_correlation", "variable": "EPAS1_vs_endothelial_score", "n": len(tumor), "groups": "tumor", "statistic": rho, "p_value": p})
    pd.DataFrame(rows).to_csv(OUT / "GSE62165_EPAS1_stage_analysis.tsv", sep="\t", index=False)
    report = ["# GSE62165 EPAS1 clinical-stage support", "", f"- Platform: GPL13667 HG-U219; samples in matrix: {len(meta)}.", f"- Mapped genes: {gene_expr.shape[0]}; endothelial program genes used: {len(up)+len(down)}.", f"- Valid grouped-stage samples: {len(valid_stage)}; available groups: {', '.join(sorted(valid_stage['grouped_stage'].unique())) if len(valid_stage) else 'none'}.", "- This cohort contains stage/tissue labels but no survival time or censoring endpoint in the downloaded series matrix.", "- Stage tests are supportive clinical association analyses only; they do not establish prognosis or causality."]
    (OUT / "GSE62165_EPAS1_stage_analysis_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("GSE62165 stage analysis complete")


if __name__ == "__main__":
    main()
