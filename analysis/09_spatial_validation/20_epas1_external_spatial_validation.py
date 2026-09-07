from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import mannwhitneyu, spearmanr


BASE = Path(r"D:\PDAC_P1")
ROOT = BASE / "analysis" / "09_spatial_validation"
RAW111 = BASE / "data" / "03_spatial" / "GSE111672" / "raw"
RAW327 = BASE / "data" / "03_spatial" / "GSE327056" / "raw"
OUT = ROOT / "EPAS1_external_validation"
OUT.mkdir(parents=True, exist_ok=True)


def read_program_genes():
    path = BASE / "analysis" / "07_candidate_screen" / "candidate_genes_stable_high_confidence_all.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df[(df["candidate_status"] == "screening_candidate") & (df["major_lineage"] == "endothelial")]
    return {
        "up": df.loc[df["direction"] == "tumor_up", "gene_id"].tolist(),
        "down": df.loc[df["direction"] == "tumor_down", "gene_id"].tolist(),
        "all": df["gene_id"].tolist(),
    }


def read_gse111(sample):
    path = RAW111 / f"GSE111672_PDAC-{sample}-indrop-filtered-expMat.txt.gz"
    df = pd.read_csv(path, sep="\t", compression="gzip")
    genes = df.iloc[:, 0].astype(str).tolist()
    labels = df.columns[1:].astype(str).tolist()
    counts = df.iloc[:, 1:].to_numpy(dtype=float)
    lib = counts.sum(axis=0)
    logcpm = np.log1p(counts / np.maximum(lib, 1.0) * 1e6)
    return genes, labels, logcpm


def single_cell_validation(program):
    rows = []
    for sample in ["A", "B"]:
        genes, labels, expr = read_gse111(sample)
        gene_index = {g: i for i, g in enumerate(genes)}
        epas = expr[gene_index["EPAS1"], :] if "EPAS1" in gene_index else np.full(expr.shape[1], np.nan)
        available = [g for g in program["all"] if g in gene_index]
        up = [gene_index[g] for g in program["up"] if g in gene_index]
        down = [gene_index[g] for g in program["down"] if g in gene_index]
        score = np.nanmean(expr[up, :], axis=0) - np.nanmean(expr[down, :], axis=0)
        labels_arr = np.asarray(labels)
        is_endo = labels_arr == "Endothelial cells"
        epas_target = epas[is_endo]
        epas_other = epas[~is_endo]
        score_target = score[is_endo]
        score_other = score[~is_endo]
        p_epas = mannwhitneyu(epas_target, epas_other, alternative="two-sided").pvalue if len(epas_target) and len(epas_other) else np.nan
        p_score = mannwhitneyu(score_target, score_other, alternative="two-sided").pvalue if len(score_target) and len(score_other) else np.nan
        rows.append({
            "dataset": "GSE111672", "sample": sample, "n_cells": len(labels_arr), "n_endothelial": int(is_endo.sum()),
            "epas1_median_endothelial": float(np.median(epas_target)) if len(epas_target) else np.nan,
            "epas1_median_other": float(np.median(epas_other)) if len(epas_other) else np.nan,
            "epas1_delta": float(np.median(epas_target) - np.median(epas_other)) if len(epas_target) else np.nan,
            "epas1_detection_endothelial": float(np.mean(epas_target > 0)) if len(epas_target) else np.nan,
            "epas1_detection_other": float(np.mean(epas_other > 0)) if len(epas_other) else np.nan,
            "epas1_wilcox_p": p_epas,
            "program_genes_used": len(available), "program_median_endothelial": float(np.median(score_target)) if len(score_target) else np.nan,
            "program_median_other": float(np.median(score_other)) if len(score_other) else np.nan,
            "program_delta": float(np.median(score_target) - np.median(score_other)) if len(score_target) else np.nan,
            "program_wilcox_p": p_score,
        })
    return pd.DataFrame(rows)


def read_visium(tag, gsm):
    prefix = f"{gsm}_{tag}"
    matrix = mmread(gzip.open(RAW327 / f"{prefix}_matrix.mtx.gz", "rb")).tocsr()
    features = pd.read_csv(gzip.open(RAW327 / f"{prefix}_features.tsv.gz", "rt"), sep="\t", header=None)
    barcodes = pd.read_csv(gzip.open(RAW327 / f"{prefix}_barcodes.tsv.gz", "rt"), sep="\t", header=None)[0].astype(str).tolist()
    positions = pd.read_csv(gzip.open(RAW327 / f"{prefix}_tissue_positions.csv.gz", "rt"))
    positions.columns = [str(c).strip() for c in positions.columns]
    if "barcode" not in positions.columns:
        positions = pd.read_csv(gzip.open(RAW327 / f"{prefix}_tissue_positions.csv.gz", "rt"), header=None)
        positions.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    positions = positions[positions["in_tissue"].astype(int) == 1]
    gene_names = features.iloc[:, 1].astype(str).tolist() if features.shape[1] > 1 else features.iloc[:, 0].astype(str).tolist()
    gene_names = pd.Index(gene_names)
    keep = [i for i, b in enumerate(barcodes) if b in set(positions["barcode"])]
    matrix = matrix[:, keep]
    kept_barcodes = [barcodes[i] for i in keep]
    positions = positions.set_index("barcode").loc[kept_barcodes].reset_index()
    lib = np.asarray(matrix.sum(axis=0)).ravel()
    epas_idx = gene_names.get_loc("EPAS1") if "EPAS1" in gene_names else None
    epas = np.log1p(np.asarray(matrix[epas_idx, :].todense()).ravel() / np.maximum(lib, 1) * 1e6) if epas_idx is not None else np.full(len(keep), np.nan)
    return positions, epas


def spatial_validation(program):
    score_path = ROOT / "GSE327056_all_spot_program_scores_global.tsv"
    scores = pd.read_csv(score_path, sep="\t")
    scaling = pd.read_csv(ROOT / "GSE327056_program_gene_global_scaling.tsv", sep="\t").set_index("gene")
    down_n = max(1, sum(g in set(program["down"]) for g in program["all"]))
    sample_info = {"A1": ("GSM9647219", "adjacent_tumor"), "B1": ("GSM9647220", "tumor"), "C1": ("GSM9647221", "tumor_stroma"), "D1": ("GSM9647222", "normal_pancreas")}
    rows = []
    for tag, (gsm, context) in sample_info.items():
        pos, epas = read_visium(tag, gsm)
        one = scores[scores["tag"] == tag].copy().set_index("barcode").loc[pos["barcode"]].reset_index()
        one["EPAS1_logCPM"] = epas
        epas_z = (epas - float(scaling.loc["EPAS1", "global_mean"])) / float(scaling.loc["EPAS1", "global_sd"])
        one["endothelial_score_exEPAS1"] = one["endothelial_score"] + epas_z / down_n
        rho, p = spearmanr(one["EPAS1_logCPM"], one["endothelial_score_exEPAS1"], nan_policy="omit")
        rows.append({"dataset": "GSE327056", "tag": tag, "context": context, "n_spots": len(one), "EPAS1_median": float(np.nanmedian(epas)), "EPAS1_detection": float(np.mean(epas > 0)), "endothelial_score_median": float(one["endothelial_score_exEPAS1"].median()), "EPAS1_endothelial_spearman_rho": float(rho), "EPAS1_endothelial_spearman_p": float(p)})
    return pd.DataFrame(rows)


def main():
    program = read_program_genes()
    sc = single_cell_validation(program)
    spatial = spatial_validation(program)
    sc.to_csv(OUT / "GSE111672_EPAS1_gene_level_validation.tsv", sep="\t", index=False)
    spatial.to_csv(OUT / "GSE327056_EPAS1_spatial_correlation.tsv", sep="\t", index=False)
    report = ["# EPAS1 external and spatial validation", "", "## GSE111672 single-cell", ""]
    for _, r in sc.iterrows():
        report.append(f"- Sample {r['sample']}: EPAS1 median {r['epas1_median_endothelial']:.3f} in endothelial vs {r['epas1_median_other']:.3f} in other cells; detection {r['epas1_detection_endothelial']:.1%} vs {r['epas1_detection_other']:.1%}; Wilcoxon P={r['epas1_wilcox_p']:.3g}.")
    report += ["", "## GSE327056 Visium", "", "- EPAS1 spot-level expression was correlated with an endothelial program score recalculated after excluding EPAS1, preventing circularity; these are localization/support analyses, not independent biological replicates.", ""]
    for _, r in spatial.iterrows():
        p_text = "<1×10⁻³⁰⁰" if r["EPAS1_endothelial_spearman_p"] == 0 else f"{r['EPAS1_endothelial_spearman_p']:.3g}"
        report.append(f"- {r['tag']} ({r['context']}): n={int(r['n_spots'])}, EPAS1 detection={r['EPAS1_detection']:.1%}, Spearman rho={r['EPAS1_endothelial_spearman_rho']:.3f}, P={p_text}.")
    report += ["", "Interpretation: EPAS1 enrichment in endothelial-labelled cells is the primary external localization result. Spatial correlations are reported with EPAS1 excluded from the program score and should be interpreted as localization support, not cell-specific causality or direct drug efficacy."]
    (OUT / "EPAS1_external_spatial_validation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("EPAS1 external/spatial validation complete")


if __name__ == "__main__":
    main()
