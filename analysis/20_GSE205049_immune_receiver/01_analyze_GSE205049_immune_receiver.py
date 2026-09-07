"""Patient-paired immune receiver sensitivity analysis for GSE205049.

This CD45-enriched paired PDAC cohort contains no endothelial compartment. It is
therefore used only to examine whether the ITGA5/ITGB1 receiver side and broad
immune-state programs differ between matched tumour and adjacent tissue.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon


ROOT = Path("D:/PDAC_P1")
DATA = ROOT / "data_external" / "GSE205049"
OUT = ROOT / "analysis" / "20_GSE205049_immune_receiver"
ANNOTATION = DATA / "GSE205049_scRNA-seq-annotations.csv.gz"
MATRIX = DATA / "GSE205049_scRNA-seq-integrated_GEM.csv.gz"
MIN_CELLS_PER_PATIENT_GROUP = 20

MYELOID_LABELS = {"Macrophages", "Monocytes", "Dendritic Cells"}
TNK_LABELS = {"CD8 T", "CD4 T", "Regulatory T", "NKT", "NK", "Proliferating T"}
RECEPTOR_GENES = ["ITGA5", "ITGB1"]
MYELOID_STATE_GENES = ["C1QA", "C1QB", "APOE", "LILRB1", "CD163", "MARCO"]
TNK_STATE_GENES = ["NKG7", "PRF1", "GZMB", "GNLY", "IFNG"]


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    """Return BH adjusted q values while preserving missing values."""
    values = p_values.to_numpy(dtype=float)
    result = np.full(values.size, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return result
    ordered = np.argsort(values[valid])
    ranked = values[valid][ordered]
    adjusted = ranked * ranked.size / np.arange(1, ranked.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(ranked)
    restored[ordered] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def read_selected_expression(matrix_path: Path, genes: list[str]) -> pd.DataFrame:
    """Stream a gene-by-cell CSV and retain only pre-specified genes."""
    selected: dict[str, np.ndarray] = {}
    gene_set = set(genes)
    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as handle:
        cell_ids = [cell_id.strip('"') for cell_id in handle.readline().rstrip("\r\n").split(",")[1:]]
        for line in handle:
            gene, separator, values = line.partition(",")
            gene = gene.strip('"')
            if separator and gene in gene_set:
                vector = np.fromstring(values, sep=",")
                if vector.size != len(cell_ids):
                    raise ValueError(f"{gene}: expected {len(cell_ids)} cells, found {vector.size}.")
                selected[gene] = vector
    missing = sorted(gene_set.difference(selected))
    if missing:
        raise ValueError(f"Pre-specified genes absent from matrix: {', '.join(missing)}")
    return pd.DataFrame(selected, index=cell_ids)


def paired_statistics(pseudobulk: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for (receiver_group, endpoint), subset in pseudobulk.groupby(["receiver_group", "endpoint"], sort=False):
        wide = subset.pivot(index="patient", columns="DiseaseState", values="value")
        if not {"PDAC", "AdjNorm"}.issubset(wide.columns):
            continue
        wide = wide.dropna(subset=["PDAC", "AdjNorm"]).copy()
        delta = wide["PDAC"] - wide["AdjNorm"]
        nonzero = delta[delta.ne(0)]
        if len(nonzero) >= 3:
            statistic, wilcoxon_p = wilcoxon(nonzero, alternative="two-sided", zero_method="wilcox", method="auto")
            sign_p = binomtest(int((nonzero > 0).sum()), n=int(len(nonzero)), p=0.5, alternative="two-sided").pvalue
        else:
            statistic, wilcoxon_p, sign_p = np.nan, np.nan, np.nan
        records.append(
            {
                "receiver_group": receiver_group,
                "endpoint": endpoint,
                "n_paired_patients": int(len(wide)),
                "n_nonzero_deltas": int(len(nonzero)),
                "median_PDAC_minus_AdjNorm": float(delta.median()),
                "mean_PDAC_minus_AdjNorm": float(delta.mean()),
                "wilcoxon_statistic": float(statistic) if np.isfinite(statistic) else np.nan,
                "paired_wilcoxon_two_sided_p": float(wilcoxon_p) if np.isfinite(wilcoxon_p) else np.nan,
                "exact_sign_test_two_sided_p": float(sign_p) if np.isfinite(sign_p) else np.nan,
                "n_positive_deltas": int((delta > 0).sum()),
                "n_negative_deltas": int((delta < 0).sum()),
            }
        )
    result = pd.DataFrame(records)
    result["BH_q"] = benjamini_hochberg(result["paired_wilcoxon_two_sided_p"])
    result["interpretation"] = np.where(
        result["BH_q"].lt(0.05),
        "paired immune-receiver sensitivity support",
        "not significant after pre-specified multiple-testing correction",
    )
    return result


def plot_paired_panels(pseudobulk: pd.DataFrame, stats: pd.DataFrame) -> None:
    endpoint_order = ["ITGA5", "ITGB1", "ITGA5/ITGB1 receptor score", "Lineage state score"]
    group_titles = {"Myeloid": "Myeloid receiver compartment", "T/NK": "T/NK receiver compartment"}
    colors = {"Myeloid": "#177E89", "T/NK": "#5B74A8"}
    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "figure.dpi": 300,
        "savefig.dpi": 450,
    })
    figure, axes = plt.subplots(2, 4, figsize=(13.5, 6.8), sharex=False)
    for row, receiver_group in enumerate(["Myeloid", "T/NK"]):
        for column, endpoint in enumerate(endpoint_order):
            axis = axes[row, column]
            subset = pseudobulk.loc[
                (pseudobulk.receiver_group == receiver_group) & (pseudobulk.endpoint == endpoint)
            ].pivot(index="patient", columns="DiseaseState", values="value")
            subset = subset.reindex(columns=["AdjNorm", "PDAC"]).dropna()
            for _, values in subset.iterrows():
                axis.plot([0, 1], values.to_numpy(), color="#B8B8B8", linewidth=0.9, zorder=1)
            axis.scatter(np.zeros(len(subset)), subset["AdjNorm"], color="#B8B8B8", s=28, zorder=2)
            axis.scatter(np.ones(len(subset)), subset["PDAC"], color=colors[receiver_group], s=32, zorder=3)
            test = stats.loc[(stats.receiver_group == receiver_group) & (stats.endpoint == endpoint)].iloc[0]
            p_value = test.paired_wilcoxon_two_sided_p
            q_value = test.BH_q
            axis.set_title(endpoint, pad=7)
            axis.set_xticks([0, 1], ["Adjacent", "PDAC"])
            axis.text(
                0.02,
                0.95,
                f"n={int(test.n_paired_patients)}\nP={p_value:.3f}; q={q_value:.3f}",
                transform=axis.transAxes,
                va="top",
                fontsize=8,
                color="#333333",
            )
            axis.spines[["top", "right"]].set_visible(False)
            axis.grid(axis="y", color="#E6E6E6", linewidth=0.7, zorder=0)
            if column == 0:
                axis.set_ylabel("Mean integrated expression")
                axis.annotate(
                    group_titles[receiver_group],
                    xy=(-0.42, 0.5),
                    xycoords="axes fraction",
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=11,
                    fontweight="bold",
                    color=colors[receiver_group],
                )
    figure.suptitle(
        "GSE205049 paired immune receiver sensitivity analysis",
        fontsize=15,
        fontweight="bold",
        y=0.985,
    )
    figure.text(
        0.5,
        0.012,
        "Patient-level pseudobulk means; only patient × tissue × lineage groups with ≥20 cells were retained. "
        "This CD45-enriched cohort does not test endothelial source expression.",
        ha="center",
        fontsize=8.5,
        color="#404040",
    )
    figure.subplots_adjust(left=0.10, right=0.99, bottom=0.12, top=0.84, wspace=0.31, hspace=0.45)
    for extension in ["png", "pdf"]:
        figure.savefig(OUT / f"Supplementary_Figure_GSE205049_immune_receiver_sensitivity.{extension}", bbox_inches="tight")
    plt.close(figure)


def write_report(annotation: pd.DataFrame, genes: list[str], pseudobulk: pd.DataFrame, stats: pd.DataFrame) -> None:
    paired = annotation.groupby(["patient", "DiseaseState"]).size().unstack(fill_value=0)
    all_pairs = int(((paired.get("PDAC", 0) > 0) & (paired.get("AdjNorm", 0) > 0)).sum())
    significant = stats.loc[stats.BH_q.lt(0.05), ["receiver_group", "endpoint", "BH_q"]]
    lines = [
        "# GSE205049 paired immune receiver sensitivity analysis",
        "",
        "## Scope and guardrail",
        "",
        "GSE205049 is a CD45-enriched, untreated, paired PDAC/adjacent-tissue single-cell cohort. "
        "It is used only to assess the immune receiver side of the ANGPT2–ITGA5/ITGB1 hypothesis. "
        "It contains no endothelial compartment and cannot validate EPAS1 endothelial localization, ANGPT2 source expression, ligand–receptor binding, or treatment response.",
        "",
        "## Design",
        "",
        f"- {all_pairs}/9 patients contained both PDAC and adjacent tissue cells.",
        f"- Immune groups were defined from the author-provided `major` labels: myeloid = macrophages/monocytes/dendritic cells; T/NK = CD4/CD8/regulatory T/NKT/NK/proliferating T.",
        f"- A patient × tissue × group required at least {MIN_CELLS_PER_PATIENT_GROUP} cells. Cells were not treated as independent replicates.",
        "- Pre-specified endpoints: ITGA5, ITGB1, their mean receptor score, and a lineage state score (myeloid macrophage program or T/NK cytotoxic program).",
        "- Two-sided paired Wilcoxon tests were used on patient pseudobulk means; Benjamini–Hochberg correction was applied across all eight tests. Exact two-sided sign tests are reported as robustness checks.",
        "",
        "## Result",
        "",
        f"- {len(significant)}/8 endpoints met BH q<0.05.",
    ]
    if len(significant):
        for row in significant.itertuples(index=False):
            lines.append(f"- {row.receiver_group} {row.endpoint}: q={row.BH_q:.3g}.")
    else:
        lines.append("- No endpoint met the pre-specified BH threshold. This cohort therefore provides no positive, independent paired support for a tumour-associated immune receiver shift; it remains a transparent sensitivity analysis rather than a negative test of the primary endothelial mechanism.")
    lines.extend([
        "",
        "## Files",
        "",
        "- `GSE205049_immune_receiver_patient_pseudobulk.tsv`: analysis-ready patient-level means and cell counts.",
        "- `GSE205049_immune_receiver_paired_statistics.tsv`: paired effect estimates, P values, sign-test P values, and BH q values.",
        "- `GSE205049_input_audit.tsv`: input alignment and retained-gene audit.",
        "- `Supplementary_Figure_GSE205049_immune_receiver_sensitivity.png/pdf`: paired-patient visualization.",
    ])
    (OUT / "GSE205049_immune_receiver_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    annotation = pd.read_csv(ANNOTATION).rename(columns={"Unnamed: 0": "cell_id", "orig.ident": "patient"})
    if annotation.cell_id.duplicated().any():
        raise ValueError("Annotation cell IDs are not unique.")
    genes = RECEPTOR_GENES + MYELOID_STATE_GENES + TNK_STATE_GENES
    expression = read_selected_expression(MATRIX, genes)
    annotation = annotation.set_index("cell_id")
    missing_annotation = expression.index.difference(annotation.index)
    missing_expression = annotation.index.difference(expression.index)
    if len(missing_annotation) or len(missing_expression):
        raise ValueError(f"Cell ID mismatch: matrix-only={len(missing_annotation)}, annotation-only={len(missing_expression)}")
    annotation = annotation.loc[expression.index].copy()
    annotation["receiver_group"] = np.select(
        [annotation.major.isin(MYELOID_LABELS), annotation.major.isin(TNK_LABELS)],
        ["Myeloid", "T/NK"],
        default="Excluded",
    )
    expression["ITGA5/ITGB1 receptor score"] = expression[RECEPTOR_GENES].mean(axis=1)
    expression["Myeloid state score"] = expression[MYELOID_STATE_GENES].mean(axis=1)
    expression["T/NK state score"] = expression[TNK_STATE_GENES].mean(axis=1)

    endpoint_map = {
        "Myeloid": ["ITGA5", "ITGB1", "ITGA5/ITGB1 receptor score", "Myeloid state score"],
        "T/NK": ["ITGA5", "ITGB1", "ITGA5/ITGB1 receptor score", "T/NK state score"],
    }
    records: list[pd.DataFrame] = []
    for receiver_group, endpoints in endpoint_map.items():
        cells = annotation.receiver_group.eq(receiver_group)
        metadata = annotation.loc[cells, ["patient", "DiseaseState"]]
        values = expression.loc[cells, endpoints].copy()
        values = values.join(metadata)
        for endpoint in endpoints:
            summary = values.groupby(["patient", "DiseaseState"], observed=True)[endpoint].agg(["mean", "size"]).reset_index()
            summary = summary.rename(columns={"mean": "value", "size": "n_cells"})
            summary["receiver_group"] = receiver_group
            summary["endpoint"] = "Lineage state score" if endpoint.endswith("state score") else endpoint
            summary = summary.loc[summary.n_cells >= MIN_CELLS_PER_PATIENT_GROUP]
            records.append(summary[["patient", "DiseaseState", "receiver_group", "endpoint", "n_cells", "value"]])
    pseudobulk = pd.concat(records, ignore_index=True)
    stats = paired_statistics(pseudobulk)
    tissue_by_patient = annotation.groupby(["patient", "DiseaseState"]).size().unstack(fill_value=0)
    paired_patient_count = int(((tissue_by_patient["PDAC"] > 0) & (tissue_by_patient["AdjNorm"] > 0)).sum())
    input_audit = pd.DataFrame([
        {"metric": "annotation_cells", "value": len(annotation)},
        {"metric": "matrix_cells", "value": len(expression)},
        {"metric": "exact_cell_id_overlap", "value": int(expression.index.equals(annotation.index))},
        {"metric": "paired_patients_with_both_tissues", "value": paired_patient_count},
        {"metric": "minimum_cells_per_patient_tissue_group", "value": MIN_CELLS_PER_PATIENT_GROUP},
        {"metric": "retained_gene_count", "value": len(genes)},
        {"metric": "retained_genes", "value": ";".join(genes)},
    ])
    pseudobulk.to_csv(OUT / "GSE205049_immune_receiver_patient_pseudobulk.tsv", sep="\t", index=False)
    stats.to_csv(OUT / "GSE205049_immune_receiver_paired_statistics.tsv", sep="\t", index=False)
    input_audit.to_csv(OUT / "GSE205049_input_audit.tsv", sep="\t", index=False)
    plot_paired_panels(pseudobulk, stats)
    write_report(annotation, genes, pseudobulk, stats)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
