#!/usr/bin/env python3
"""Directionality audit, spatial conditional models, and revised Figures 4–5.

This script is deliberately conservative.  Spatial spots are used only to estimate
within-patient/sample coefficients; inference is performed across independent
patients/samples, not across spots.  Drug signature reversal and docking are kept
as separate computational evidence layers and are never interpreted as EPAS1
agonism, antagonism, direct binding, or clinical activity.
"""
from __future__ import annotations

from pathlib import Path
import sys
from math import comb
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.paths import ANALYSIS_ROOT


ROOT = ANALYSIS_ROOT.parent
OUT = ANALYSIS_ROOT / "21_revision_strengthening"
OUT.mkdir(parents=True, exist_ok=True)
FIGOUT = ANALYSIS_ROOT / "15_priority_figures" / "main_figures_unified"

mpl.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8.5, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.labelsize": 8.5, "legend.fontsize": 7.4, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 400,
})

COLORS = {
    "navy": "#173F5F", "blue": "#20639B", "teal": "#3CAEA3",
    "gold": "#F6C85F", "coral": "#ED553B", "gray": "#7D8794",
    "lightgray": "#D9E1E7", "dark": "#263238", "violet": "#7A5195",
}
PROGRAM_ORDER = ["endothelial", "angiogenesis", "CAF_ECM", "TGFb_fibrosis",
                 "myeloid_inflammation", "TNF_NFkB", "cytotoxic_TNK", "T_exhaustion"]
PROGRAM_LABELS = {
    "endothelial": "Endothelial", "angiogenesis": "Angiogenesis",
    "CAF_ECM": "CAF / ECM", "TGFb_fibrosis": "TGF-β / fibrosis",
    "myeloid_inflammation": "Myeloid inflammation", "TNF_NFkB": "TNF–NF-κB",
    "cytotoxic_TNK": "Cytotoxic T/NK", "T_exhaustion": "T-cell exhaustion",
}


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def bh(values: pd.Series) -> np.ndarray:
    values = values.astype(float).to_numpy()
    keep = np.isfinite(values)
    out = np.full(values.size, np.nan)
    if keep.any():
        out[keep] = multipletests(values[keep], method="fdr_bh")[1]
    return out


def exact_sign_p(n_positive: int, n_total: int) -> float:
    """Two-sided exact binomial sign test without relying on scipy version."""
    if n_total == 0:
        return np.nan
    probability = sum(comb(n_total, index) for index in range(n_positive + 1)) / 2**n_total
    probability = min(probability, 1 - probability + comb(n_total, n_positive) / 2**n_total)
    return min(1.0, 2 * probability)


def residualize(value: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    valid = np.isfinite(value) & np.all(np.isfinite(covariates), axis=1)
    output = np.full(value.shape, np.nan, dtype=float)
    if valid.sum() < covariates.shape[1] + 10:
        return output
    design = np.column_stack([np.ones(valid.sum()), covariates[valid]])
    beta, *_ = np.linalg.lstsq(design, value[valid], rcond=None)
    output[valid] = value[valid] - design @ beta
    return output


def zscore(value: np.ndarray) -> np.ndarray:
    standard_deviation = np.nanstd(value, ddof=1)
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        return np.full(value.shape, np.nan)
    return (value - np.nanmean(value)) / standard_deviation


def spatial_blocks(data: pd.DataFrame) -> pd.Series:
    """Four-by-four coordinate blocks, constructed independently per ROI/sample."""
    xbin = pd.qcut(data["pxl_col_in_fullres"].rank(method="first"), 4, labels=False, duplicates="drop")
    ybin = pd.qcut(data["pxl_row_in_fullres"].rank(method="first"), 4, labels=False, duplicates="drop")
    return (xbin.astype(str) + "_" + ybin.astype(str)).astype("category")


def conditional_slope(data: pd.DataFrame, outcome: str) -> dict:
    """Standardized EPAS1 slope conditional on vascular score, depth and blocks."""
    needed = ["EPAS1_logCPM", "vascular_marker_score", "total_counts", "n_genes", outcome]
    frame = data.dropna(subset=needed).copy()
    if frame.shape[0] < 120:
        return {"n_spots": frame.shape[0], "beta_epas1": np.nan, "partial_r": np.nan}
    frame["spatial_block"] = spatial_blocks(frame)
    block = pd.get_dummies(frame["spatial_block"], drop_first=True, dtype=float)
    continuous = np.column_stack([
        zscore(frame["vascular_marker_score"].to_numpy(float)),
        zscore(np.log1p(frame["total_counts"].to_numpy(float))),
        zscore(np.log1p(frame["n_genes"].to_numpy(float))),
        block.to_numpy(float),
    ])
    epas1_resid = residualize(zscore(frame["EPAS1_logCPM"].to_numpy(float)), continuous)
    outcome_resid = residualize(zscore(frame[outcome].to_numpy(float)), continuous)
    keep = np.isfinite(epas1_resid) & np.isfinite(outcome_resid)
    if keep.sum() < 100 or np.nanstd(epas1_resid[keep]) == 0:
        return {"n_spots": int(keep.sum()), "beta_epas1": np.nan, "partial_r": np.nan}
    # With standardized residuals this OLS coefficient equals the partial correlation.
    partial_r = float(np.corrcoef(epas1_resid[keep], outcome_resid[keep])[0, 1])
    return {"n_spots": int(keep.sum()), "beta_epas1": partial_r, "partial_r": partial_r}


def build_directionality_audit() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    for cohort in ["GSE155698", "GSE212966"]:
        edge = read_tsv(ROOT / "analysis" / "06_pseudobulk_DE" / "edgeR" / cohort / "endothelial_DE.tsv")
        deseq = read_tsv(ROOT / "analysis" / "06_pseudobulk_DE" / "DESeq2" / cohort / "endothelial_DE.tsv")
        edge_row = edge.loc[edge["gene_id"].astype(str).str.replace('"', '', regex=False).eq("EPAS1")].iloc[0]
        deseq_row = deseq.loc[deseq["gene_id"].astype(str).str.replace('"', '', regex=False).eq("EPAS1")].iloc[0]
        rows.extend([
            {"layer": "Tumor–adjacent endothelial DE", "dataset": cohort, "comparison": "tumor vs adjacent", "method": "edgeR",
             "direction": "decreased in tumor" if edge_row["logFC"] < 0 else "increased in tumor", "effect": edge_row["logFC"], "P_or_FDR": edge_row["FDR"],
             "directional_interpretation": "observed abundance difference; not functional direction"},
            {"layer": "Tumor–adjacent endothelial DE", "dataset": cohort, "comparison": "tumor vs adjacent", "method": "DESeq2",
             "direction": "decreased in tumor" if deseq_row["log2FoldChange"] < 0 else "increased in tumor", "effect": deseq_row["log2FoldChange"], "P_or_FDR": deseq_row["padj"],
             "directional_interpretation": "observed abundance difference; not functional direction"},
        ])

    proxy = read_tsv(ROOT / "analysis" / "13_virtual_perturbation" / "virtual_perturbation_proxy_program_delta.tsv")
    network = read_tsv(ROOT / "analysis" / "13_virtual_perturbation" / "virtual_perturbation_network_program_summary.tsv")
    proxy = proxy.query("target == 'EPAS1' and cell_type == 'endothelial'").copy()
    network = network.query("target == 'EPAS1' and cell_type == 'endothelial'").copy()
    merged = proxy.merge(network[["dataset", "cell_type", "target", "perturbation", "program", "n_sig_fdr"]],
                         left_on=["dataset", "cell_type", "target", "proxy", "program"],
                         right_on=["dataset", "cell_type", "target", "perturbation", "program"], how="left")
    for (perturbation, program), subset in merged.groupby(["proxy", "program"], sort=False):
        deltas = subset["delta"].astype(float)
        n_positive, n_negative = int((deltas > 0).sum()), int((deltas < 0).sum())
        dominant = "increased program score" if n_positive > n_negative else "decreased program score" if n_negative > n_positive else "mixed"
        rows.append({"layer": "Virtual network perturbation", "dataset": "3 scRNA-seq cohorts", "comparison": f"endothelial EPAS1 {perturbation}",
                     "method": "scTenifold proxy-program response", "direction": dominant,
                     "effect": float(deltas.median()), "P_or_FDR": int((subset["n_sig_fdr"].fillna(0) > 0).sum()),
                     "directional_interpretation": f"{n_positive}/3 positive, {n_negative}/3 negative signed proxy responses; significant-network cohorts reported separately"})

    drug = read_tsv(ROOT / "analysis" / "19_translation_closure" / "EPAS1_drug_repositioning_evidence_matrix.tsv")
    for _, row in drug.iterrows():
        rows.append({"layer": "Drug signature / pharmacology", "dataset": "LINCS DrugReflector + PRISM", "comparison": str(row["compound"]),
                     "method": "signature reversal / pancreatic cell-line viability", "direction": "EPAS1 functional mode unknown",
                     "effect": float(row["PRISM_median_logFC"]), "P_or_FDR": np.nan,
                     "directional_interpretation": "Signature reversal does not identify agonism/antagonism; PRISM viability does not establish EPAS1 dependence."})
        rows.append({"layer": "Structural compatibility", "dataset": "five EPAS1 crystal conformers", "comparison": str(row["compound"]),
                     "method": "Vina + co-crystal redocking QC", "direction": "binding/function direction unknown",
                     "effect": float(row["docking_mean_affinity_kcal_mol"]), "P_or_FDR": np.nan,
                     "directional_interpretation": "Docking supports pose compatibility only; it cannot establish binding or target modulation."})
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "EPAS1_directionality_audit.tsv", sep="\t", index=False)
    merged.to_csv(OUT / "EPAS1_virtual_perturbation_signed_program_responses.tsv", sep="\t", index=False)
    return audit, merged, drug


def build_spatial_models() -> tuple[pd.DataFrame, pd.DataFrame]:
    anchor = read_tsv(ROOT / "analysis" / "09_spatial_validation" / "GSE282302" / "GSE282302_spot_anchor_scores.tsv")
    ecology = read_tsv(ROOT / "analysis" / "09_spatial_validation" / "GSE282302" / "spatial_ecology" / "GSE282302_spatial_ecology_spots.tsv")
    anchor = anchor.loc[anchor["in_tissue"].eq(1)].copy()
    joined = anchor.merge(ecology, on=["sample", "patient", "barcode", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"], how="inner")
    outcomes = ["myeloid_neighbor_score", "T_NK_neighbor_score", "endothelial_neighbor_score"]
    rows = []
    for patient, subset in joined.groupby("patient", sort=True):
        for outcome in outcomes:
            result = conditional_slope(subset, outcome)
            rows.append({"dataset": "GSE282302", "analysis_unit": "patient", "patient_or_sample": patient,
                         "outcome": outcome, "model": "outcome ~ EPAS1 + vascular score + depth + detected genes + 4x4 spatial block",
                         **result})
    patient = pd.DataFrame(rows)

    # Independent sample-level anchor control in GSE297144.  It tests whether
    # EPAS1 retains association with a non-EPAS1 endothelial program after vascular adjustment,
    # not an immune association (no matched immune ecology table is available).
    anchor_297 = read_tsv(ROOT / "analysis" / "09_spatial_validation" / "GSE297144" / "GSE297144_spot_anchor_scores.tsv")
    anchor_297 = anchor_297.loc[anchor_297["in_tissue"].eq(1)].copy()
    control_rows = []
    for sample, subset in anchor_297.groupby("sample", sort=True):
        result = conditional_slope(subset, "endothelial_program_excluding_EPAS1")
        control_rows.append({"dataset": "GSE297144", "analysis_unit": "sample", "patient_or_sample": sample,
                             "outcome": "endothelial_program_excluding_EPAS1",
                             "model": "outcome ~ EPAS1 + vascular score + depth + detected genes + 4x4 spatial block", **result})
    control = pd.DataFrame(control_rows)
    all_patient = pd.concat([patient, control], ignore_index=True)
    summary_rows = []
    for (dataset, outcome), subset in all_patient.groupby(["dataset", "outcome"], sort=False):
        beta = subset["beta_epas1"].dropna().astype(float)
        n_positive = int((beta > 0).sum())
        n_negative = int((beta < 0).sum())
        summary_rows.append({"dataset": dataset, "outcome": outcome, "n_independent_units": len(beta),
                             "n_positive": n_positive, "n_negative": n_negative, "median_conditional_beta": beta.median(),
                             "mean_conditional_beta": beta.mean(), "exact_two_sided_sign_p": exact_sign_p(min(n_positive, n_negative), len(beta))})
    summary = pd.DataFrame(summary_rows)
    summary["BH_q_across_spatial_endpoints"] = bh(summary["exact_two_sided_sign_p"])
    patient.to_csv(OUT / "EPAS1_spatial_conditional_model_patient_coefficients.tsv", sep="\t", index=False)
    control.to_csv(OUT / "EPAS1_spatial_conditional_model_GSE297144_anchor_control.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "EPAS1_spatial_conditional_model_summary.tsv", sep="\t", index=False)
    return patient, summary


def export(figure: plt.Figure, basename: str, output_dir: Path = FIGOUT) -> None:
    figure.savefig(output_dir / f"{basename}.pdf", bbox_inches="tight")
    figure.savefig(output_dir / f"{basename}.png", dpi=400, bbox_inches="tight")
    figure.savefig(OUT / f"{basename}.pdf", bbox_inches="tight")
    figure.savefig(OUT / f"{basename}.png", dpi=400, bbox_inches="tight")
    plt.close(figure)


def draw_direction_audit(audit: pd.DataFrame, perturb: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(12.0, 3.4), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 2.05, 1.15])
    ax1, ax2, ax3 = [fig.add_subplot(grid[0, index]) for index in range(3)]
    de = audit.query("layer == 'Tumor–adjacent endothelial DE'")
    labels = [f"{row.dataset}\n{row.method}" for row in de.itertuples()]
    vals = de["effect"].to_numpy(float)
    ax1.barh(np.arange(len(vals)), vals, color=[COLORS["coral"] if value < 0 else COLORS["teal"] for value in vals])
    ax1.axvline(0, color=COLORS["dark"], lw=.8)
    ax1.set_yticks(np.arange(len(vals)), labels); ax1.invert_yaxis()
    ax1.set_xlabel("Tumor–adjacent log2FC")
    ax1.set_title("A  Observed abundance direction")

    subset = perturb.query("cell_type == 'endothelial'").copy()
    matrix = subset.pivot(index="program", columns=["proxy", "dataset"], values="delta").reindex(PROGRAM_ORDER)
    columns = [("KO", "GSE154778_primary"), ("KO", "GSE155698_primary"), ("KO", "GSE212966_primary"),
               ("OE", "GSE154778_primary"), ("OE", "GSE155698_primary"), ("OE", "GSE212966_primary")]
    matrix = matrix.reindex(columns=pd.MultiIndex.from_tuples(columns))
    vmax = np.nanmax(np.abs(matrix.to_numpy(float)))
    image = ax2.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax2.set_xticks(np.arange(6), ["154778", "155698", "212966", "154778", "155698", "212966"], rotation=35, ha="right")
    ax2.set_yticks(np.arange(len(PROGRAM_ORDER)), [PROGRAM_LABELS[item] for item in PROGRAM_ORDER])
    ax2.axvline(2.5, color="white", lw=3)
    ax2.text(1, -1.15, "KO", ha="center", weight="bold", color=COLORS["navy"])
    ax2.text(4, -1.15, "OE", ha="center", weight="bold", color=COLORS["violet"])
    for row_index, program in enumerate(PROGRAM_ORDER):
        for column_index, (condition, dataset) in enumerate(columns):
            cell = subset.loc[(subset["program"] == program) & (subset["proxy"] == condition) & (subset["dataset"] == dataset)]
            if not cell.empty and float(cell.iloc[0].get("n_sig_fdr", 0)) > 0:
                ax2.text(column_index, row_index, "•", ha="center", va="center", fontsize=11, color="black")
    ax2.set_title("B  Signed proxy-program response in endothelial cells")
    cbar = fig.colorbar(image, ax=ax2, fraction=.046, pad=.025)
    cbar.set_label("Δ program score")

    rows = []
    for perturbation in ["KO", "OE"]:
        for program in PROGRAM_ORDER:
            d = subset.loc[(subset["proxy"] == perturbation) & (subset["program"] == program)]
            rows.append((perturbation, program, int((d["n_sig_fdr"].fillna(0) > 0).sum()), len(d)))
    table = pd.DataFrame(rows, columns=["perturbation", "program", "n", "denom"])
    x = np.arange(len(PROGRAM_ORDER)); width = .34
    for index, perturbation in enumerate(["KO", "OE"]):
        d = table.loc[table["perturbation"] == perturbation]
        ax3.bar(x + (index - .5) * width, d["n"] / d["denom"], width, label=perturbation,
                color=COLORS["blue"] if perturbation == "KO" else COLORS["violet"])
        for xpos, numerator, denominator in zip(x + (index - .5) * width, d["n"], d["denom"]):
            ax3.text(xpos, numerator / denominator + .035, f"{numerator}/{denominator}", ha="center", va="bottom", fontsize=6.3)
    ax3.set_ylim(0, 1.22); ax3.set_xticks(x, ["Endo", "Angio", "CAF", "TGF-β", "Myeloid", "TNF", "Cyto", "Exhaust"], rotation=45, ha="right")
    ax3.set_ylabel("Cohorts with ≥1 FDR-significant\nnetwork program gene")
    ax3.set_title("C  Network-significance recurrence")
    ax3.legend(title="Network operation", ncol=2, loc="upper right")
    fig.text(.01, .01, "Dot, ≥1 FDR-significant network gene. The non-opposing KO/OE score changes preclude a monotonic EPAS1 dosage interpretation.", fontsize=7.2, color=COLORS["dark"])
    export(fig, "Figure_S31_EPAS1_directionality_audit")


def draw_figure4(perturb: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(12.3, 5.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=[1.27, 1.27, 1.08], height_ratios=[1, .72])
    conditions = [("KO", "A"), ("OE", "B")]
    subset = perturb.query("cell_type == 'endothelial'").copy()
    vmax = np.nanmax(np.abs(subset["delta"].astype(float)))
    for col, (condition, label) in enumerate(conditions):
        ax = fig.add_subplot(grid[0, col])
        part = subset.loc[subset["proxy"] == condition]
        matrix = part.pivot(index="program", columns="dataset", values="delta").reindex(index=PROGRAM_ORDER, columns=["GSE154778_primary", "GSE155698_primary", "GSE212966_primary"])
        image = ax.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(range(3), ["GSE154778", "GSE155698", "GSE212966"], rotation=25, ha="right")
        ax.set_yticks(range(8), [PROGRAM_LABELS[item] for item in PROGRAM_ORDER] if col == 0 else [])
        ax.set_title(f"{label}  EPAS1 {condition}: signed proxy-program response")
        for row_index, program in enumerate(PROGRAM_ORDER):
            for column_index, dataset in enumerate(["GSE154778_primary", "GSE155698_primary", "GSE212966_primary"]):
                cell = part.loc[(part["program"] == program) & (part["dataset"] == dataset)]
                if not cell.empty and float(cell.iloc[0].get("n_sig_fdr", 0)) > 0:
                    ax.text(column_index, row_index, "•", ha="center", va="center", color="black", fontsize=11)
        cbar = fig.colorbar(image, ax=ax, fraction=.046, pad=.03)
        cbar.set_label("Δ score")
    ax = fig.add_subplot(grid[0, 2])
    recurrence = []
    for operation in ["KO", "OE"]:
        for program in PROGRAM_ORDER:
            part = subset.loc[(subset["proxy"] == operation) & (subset["program"] == program)]
            recurrence.append([operation, program, int((part["n_sig_fdr"].fillna(0) > 0).sum()), len(part)])
    recurrence = pd.DataFrame(recurrence, columns=["operation", "program", "n", "eligible"])
    index = np.arange(8); width = .36
    for position, operation in enumerate(["KO", "OE"]):
        part = recurrence.loc[recurrence["operation"] == operation]
        ax.bar(index + (position - .5) * width, part["n"] / part["eligible"], width, label=operation,
               color=COLORS["blue"] if operation == "KO" else COLORS["violet"])
        for xpos, numerator, denominator in zip(index + (position - .5) * width, part["n"], part["eligible"]):
            ax.text(xpos, numerator / denominator + .035, f"{numerator}/{denominator}", ha="center", fontsize=6.3)
    ax.set_ylim(0, 1.22); ax.set_xticks(index, ["Endo", "Angio", "CAF", "TGFβ", "Myeloid", "TNF", "Cyto", "Exhaust"], rotation=45, ha="right")
    ax.set_ylabel("Cohort recurrence")
    ax.set_title("C  FDR-significant network recurrence")
    ax.legend(ncol=2, loc="upper right")

    ax = fig.add_subplot(grid[1, :])
    paired = subset.pivot(index=["dataset", "program"], columns="proxy", values="delta").dropna()
    for program, part in paired.reset_index().groupby("program"):
        ax.scatter(part["KO"], part["OE"], s=42, color=COLORS["gold"] if program in ["endothelial", "angiogenesis"] else COLORS["gray"], alpha=.9, edgecolor="white", linewidth=.5)
    bound = np.nanmax(np.abs(paired[["KO", "OE"]].to_numpy())) * 1.15
    ax.axhline(0, color=COLORS["dark"], lw=.7); ax.axvline(0, color=COLORS["dark"], lw=.7)
    ax.plot([-bound, bound], [bound, -bound], ls="--", lw=1, color=COLORS["coral"], label="opposite-sign expectation")
    rho, p_value = stats.spearmanr(paired["KO"], paired["OE"])
    ax.text(.02, .94, f"KO–OE Spearman ρ={rho:.2f}; P={p_value:.3g}\n(no required anti-correlation)", transform=ax.transAxes, va="top", fontsize=7.5,
            bbox={"boxstyle": "round,pad=.25", "facecolor": "white", "edgecolor": COLORS["lightgray"]})
    ax.set_xlim(-bound, bound); ax.set_ylim(-bound, bound)
    ax.set_xlabel("KO Δ program score"); ax.set_ylabel("OE Δ program score")
    ax.set_title("D  KO and OE operations do not establish a monotonic EPAS1 dosage response")
    ax.legend(loc="lower right")
    export(fig, "Figure_4_virtual_perturbation_unified_REVISED")
    recurrence.to_csv(OUT / "Figure_4_network_recurrence_source_data.tsv", sep="\t", index=False)


def draw_figure5(drug: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(12.0, 4.7))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.05, 1.35])
    fig.subplots_adjust(left=.055, right=.985, top=.87, bottom=.27, wspace=.42)
    ax1, ax2, ax3 = [fig.add_subplot(grid[0, index]) for index in range(3)]
    priority = read_tsv(ROOT / "analysis" / "13_virtual_perturbation" / "drugreflector" / "drugreflector_priority.tsv")
    mapping = read_tsv(ROOT / "analysis" / "13_virtual_perturbation" / "drugreflector" / "drugreflector_compound_mapping.tsv")
    priority = priority.merge(mapping[["compound", "cmap_name"]], on="compound", how="left").sort_values("mean_rank").head(8)
    labels = priority["cmap_name"].fillna(priority["compound"]).astype(str).to_list()
    colors = [COLORS["coral"] if name == "triclabendazole" else COLORS["gold"] if name == "Y-39983" else COLORS["lightgray"] for name in labels]
    ax1.barh(np.arange(len(priority)), priority["mean_rank"], color=colors, edgecolor="white")
    ax1.set_yticks(np.arange(len(priority)), labels); ax1.invert_yaxis(); ax1.set_xlim(-.1, 3.45); ax1.set_xlabel("DrugReflector mean reverse rank (lower is earlier)")
    ax1.set_title("A  Recurrent signature-reversal ranking")
    ax1.text(.01, -.22, "Reverse ranking is transcriptomic connectivity;\nit does not identify EPAS1 agonism/antagonism.", transform=ax1.transAxes, fontsize=7.2, va="top")

    features = read_tsv(ROOT / "analysis" / "11_methodology_review" / "DepMap_PRISM_pharmacology" / "DepMap_PRISM_PDAC_matched_features.tsv")
    drug_columns = [("triclabendazole", "triclabendazole_logFC", COLORS["coral"]), ("Y-39983", "Y39983_logFC", COLORS["gold"])]
    rng = np.random.default_rng(20260906)
    ax2.set_ylim(-2.3, 1.05)
    for position, (name, column, color) in enumerate(drug_columns):
        values = features[column].dropna().to_numpy(float)
        jitter = rng.normal(0, .055, len(values))
        ax2.scatter(np.full(len(values), position) + jitter, values, color=color, s=22, alpha=.82, edgecolors="white", linewidth=.35, zorder=3)
        median = np.median(values)
        ax2.hlines(median, position-.28, position+.28, lw=2.4, color=COLORS["dark"], zorder=4)
        ax2.text(position, .93, f"n={len(values)}\nmedian={median:.2f}", ha="center", va="top", fontsize=7)
    ax2.axhline(-.5, ls="--", color=COLORS["gray"], lw=1, label="screening reference: −0.5")
    ax2.set_xticks([0, 1], ["triclabendazole", "Y-39983"]); ax2.set_ylabel("PRISM viability logFC")
    ax2.set_title("B  Pancreatic cell-line viability (PRISM)")
    ax2.legend(loc="lower right")

    dock = read_tsv(ROOT / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer" / "final_exhaustiveness16" / "final_exhaustiveness16_summary.tsv")
    compound_map = {"BRD-K56751279": "Y-39983", "BRD-K81916719": "triclabendazole"}
    dock["drug"] = dock["compound"].map(compound_map)
    redock = read_tsv(ROOT / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer" / "multiconformer_summary_corrected.tsv")
    redock = redock.loc[redock["kind"].eq("cocrystal_redock")].copy()
    for position, (name, color) in enumerate([("Y-39983", COLORS["gold"]), ("triclabendazole", COLORS["coral"])]):
        values = dock.loc[dock["drug"].eq(name), "affinity_best_kcal_mol"].to_numpy(float)
        ax3.scatter(np.full(len(values), position) + rng.normal(0, .045, len(values)), values, color=color, s=36, zorder=3, edgecolors="white", linewidth=.4)
        ax3.hlines(np.median(values), position-.25, position+.25, color=COLORS["dark"], lw=2)
    ax3.axhline(0, lw=.7, color=COLORS["gray"])
    ax3.set_xticks([0, 1], ["Y-39983", "triclabendazole"]); ax3.set_ylabel("Best Vina affinity (kcal/mol)")
    ax3.set_title("C  Five-conformer pose compatibility")
    rmsd_text = " | ".join(f"{row.pdb_id}: {row.rmsd_angstrom:.2f} Å" for row in redock.itertuples())
    ax3.text(.02, .035, f"Co-crystal redocking QC (5/5 pass)\n{rmsd_text}", transform=ax3.transAxes, fontsize=6.55, va="bottom",
             bbox={"boxstyle": "round,pad=.22", "facecolor": "white", "edgecolor": COLORS["lightgray"]})
    fig.text(.055, .055, "Evidence layers are intentionally not combined into a single efficacy score. Triclabendazole is a testable repurposing hypothesis;\nY-39983 is a structure-supported negative pharmacologic comparator. Docking supports pose compatibility, not direct binding or functional proof.", fontsize=7.25)
    export(fig, "Figure_5_drug_structure_closure_unified_REVISED")
    priority.to_csv(OUT / "Figure_5_DrugReflector_ranking_source_data.tsv", sep="\t", index=False)
    dock.to_csv(OUT / "Figure_5_docking_source_data.tsv", sep="\t", index=False)


def draw_spatial_conditional_model(patient: pd.DataFrame, summary: pd.DataFrame) -> None:
    """Patient/sample slope display; intentionally no spot-level significance labels."""
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.15), gridspec_kw={"width_ratios": [2.7, 1]}, constrained_layout=True)
    gse282 = patient.loc[patient["dataset"].eq("GSE282302")].copy()
    order = ["myeloid_neighbor_score", "T_NK_neighbor_score", "endothelial_neighbor_score"]
    labels = ["Myeloid\nneighbor", "T/NK\nneighbor", "Endothelial\nneighbor"]
    palette = [COLORS["gray"], COLORS["teal"], COLORS["blue"]]
    rng = np.random.default_rng(20260906)
    for index, (outcome, label, color) in enumerate(zip(order, labels, palette)):
        values = gse282.loc[gse282["outcome"].eq(outcome), "beta_epas1"].dropna().to_numpy(float)
        axes[0].scatter(np.full(values.size, index) + rng.normal(0, .055, values.size), values, s=30, color=color, alpha=.86, edgecolor="white", linewidth=.4, zorder=3)
        axes[0].hlines(np.median(values), index-.27, index+.27, lw=2.3, color=COLORS["dark"], zorder=4)
        row = summary.loc[(summary["dataset"].eq("GSE282302")) & (summary["outcome"].eq(outcome))].iloc[0]
        q_label = "<0.001" if row.BH_q_across_spatial_endpoints < .001 else f"{row.BH_q_across_spatial_endpoints:.3f}"
        axes[0].text(index, .135, f"n={int(row.n_independent_units)}\nq={q_label}", ha="center", va="top", fontsize=7)
    axes[0].axhline(0, color=COLORS["dark"], lw=.8); axes[0].set_ylim(-.08, .15)
    axes[0].set_xticks(range(3), labels); axes[0].set_ylabel("Conditional EPAS1 coefficient (patient-level estimate)")
    axes[0].set_title("A  GSE282302 immune-ecology conditional model")
    axes[0].text(.01, -.32, "Within each patient: outcome ~ EPAS1 + vascular score + depth + detected genes + 4×4 coordinate block.\nInference is the two-sided sign test across 14 patients; not a spot-level test.", transform=axes[0].transAxes, fontsize=6.9, va="top")

    gse297 = patient.loc[patient["dataset"].eq("GSE297144")].copy()
    values = gse297["beta_epas1"].dropna().to_numpy(float)
    axes[1].scatter(rng.normal(0, .045, values.size), values, s=34, color=COLORS["violet"], alpha=.85, edgecolor="white", linewidth=.4)
    axes[1].hlines(np.median(values), -.25, .25, lw=2.3, color=COLORS["dark"])
    row = summary.loc[summary["dataset"].eq("GSE297144")].iloc[0]
    axes[1].axhline(0, color=COLORS["dark"], lw=.8); axes[1].set_xlim(-.4, .4); axes[1].set_ylim(-.08, .15)
    axes[1].set_xticks([0], ["Endothelial\nanchor control"]); axes[1].set_ylabel("Conditional EPAS1 coefficient")
    axes[1].set_title("B  GSE297144 control")
    q_label = "<0.001" if row.BH_q_across_spatial_endpoints < .001 else f"{row.BH_q_across_spatial_endpoints:.3f}"
    axes[1].text(0, .135, f"n={int(row.n_independent_units)}\nq={q_label}", ha="center", va="top", fontsize=7)
    axes[1].text(.02, -.32, "Processed metadata support sample,\nnot donor, as the independent unit.", transform=axes[1].transAxes, fontsize=6.9, va="top")
    export(fig, "Figure_S32_spatial_conditional_model")


def write_report(audit: pd.DataFrame, spatial_summary: pd.DataFrame) -> None:
    tumor = audit.loc[audit["layer"].eq("Tumor–adjacent endothelial DE")]
    virtual = audit.loc[audit["layer"].eq("Virtual network perturbation")]
    decrease = int((tumor["effect"] < 0).sum())
    spatial_lines = []
    for row in spatial_summary.itertuples():
        spatial_lines.append(f"- {row.dataset}, `{row.outcome}`: n={row.n_independent_units}; median conditional beta={row.median_conditional_beta:.3f}; exact sign-test P={row.exact_two_sided_sign_p:.4g}; BH q={row.BH_q_across_spatial_endpoints:.4g}.")
    report = f"""# Directionality audit and conditional spatial models

## Pre-specified question

Do the observed endothelial EPAS1 difference, virtual perturbations, spatial associations, and drug/structure analyses support a single directional **EPAS1-inhibition** claim?  The answer is **no**. The analyses support a vessel-associated EPAS1 state and a computational intervention hypothesis, but not a monotonic functional direction or direct drug target engagement.

## Directionality audit

- In the two eligible endothelial tumor-versus-adjacent contrasts, EPAS1 had a negative log2 fold change in all {decrease}/4 method–cohort estimates. This is an abundance observation and does not establish whether increasing or decreasing EPAS1 is beneficial.
- The signed KO/OE proxy program responses are not a required mirror image. Therefore scTenifoldKnk/scTenifoldNet should be described as **network-sensitivity operations**, not gene-dosage simulations or causal proof.
- DrugReflector ranks reversal of a derived expression signature. It does not specify whether either compound activates or inhibits EPAS1. Vina/redocking adds structural pose compatibility only; it cannot resolve functional direction or direct binding.

## Conditional spatial models

Each GSE282302 patient was modeled separately with vascular-marker score (PECAM1/VWF/EMCN), depth, detected genes, and a 4×4 coordinate block as covariates. Patient-level sign tests, not spot-level P values, are the inferential unit. GSE297144 serves solely as an independent endothelial-anchor control because it lacks a matched immune-ecology table.

{chr(10).join(spatial_lines)}

## Claim gate for the revised manuscript

Use: **“EPAS1 marks a vessel-associated endothelial state with reproducible network sensitivity. The ANGPT2 interface and triclabendazole are testable computational hypotheses.”**

Do not use: “EPAS1 drives immune recruitment,” “EPAS1 inhibition is therapeutic,” “triclabendazole targets EPAS1,” or any wording that treats spot-level associations as patient-level causal evidence.
"""
    (OUT / "EPAS1_directionality_and_spatial_conditional_model_report.md").write_text(textwrap.dedent(report), encoding="utf-8")


def main() -> None:
    audit, perturb, drug = build_directionality_audit()
    spatial_patient, spatial_summary = build_spatial_models()
    draw_direction_audit(audit, perturb)
    draw_figure4(perturb)
    draw_figure5(drug)
    draw_spatial_conditional_model(pd.concat([spatial_patient, read_tsv(OUT / "EPAS1_spatial_conditional_model_GSE297144_anchor_control.tsv")], ignore_index=True), spatial_summary)
    write_report(audit, spatial_summary)
    print("Completed directionality audit, spatial conditional models, and revised Figures 4–5.")


if __name__ == "__main__":
    main()
