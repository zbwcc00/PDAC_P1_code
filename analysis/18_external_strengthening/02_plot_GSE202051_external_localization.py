from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


ROOT = Path("D:/PDAC_P1")
OUT = ROOT / "analysis" / "18_external_strengthening"
H5AD = ROOT / "data_external" / "GSE202051" / "GSE202051_totaldata-final-toshare.h5ad"

LINEAGE_COLORS = {
    "Epithelial": "#D65F5F", "Fibroblast": "#C68B59", "Immune": "#5E9ACF",
    "Endothelial": "#177E89", "Endocrine": "#8F6BB3", "Schwann": "#6FAF72", "unknown": "#A6A6A6",
}


def panel_label(axis, label: str) -> None:
    axis.text(-0.08, 1.04, label, transform=axis.transAxes, fontsize=15, fontweight="bold", va="bottom")


def main() -> None:
    np.random.seed(20260906)
    summary = pd.read_csv(OUT / "GSE202051_independent_localization_summary.tsv", sep="\t").iloc[0]
    patient_lineage = pd.read_csv(OUT / "GSE202051_patient_lineage_expression.tsv", sep="\t")
    subtype = pd.read_csv(OUT / "GSE202051_endothelial_subtype_expression.tsv", sep="\t")
    data = ad.read_h5ad(H5AD, backed="r")
    obs = data.obs
    coordinates = np.asarray(data.obsm["X_umap"])
    epas_idx = data.var_names.get_loc("EPAS1")
    epas = data[:, epas_idx].X.toarray().ravel()
    lineage = obs["broad_celltypes"].astype(str).to_numpy()

    # Preserve endothelial cells and randomly downsample all other cells for legibility.
    endo = lineage == "Endothelial"
    other_index = np.flatnonzero(~endo)
    sampled_other = np.random.choice(other_index, size=min(80000, len(other_index)), replace=False)
    display_index = np.concatenate([sampled_other, np.flatnonzero(endo)])
    display_lineage = lineage[display_index]

    figure = plt.figure(figsize=(15, 10.5), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=[1.05, 0.95], width_ratios=[1.05, 0.95])
    axis_a = figure.add_subplot(grid[0, 0])
    for group in ["Epithelial", "Fibroblast", "Immune", "Endocrine", "Schwann", "unknown", "Endothelial"]:
        mask = display_lineage == group
        if np.any(mask):
            axis_a.scatter(coordinates[display_index[mask], 0], coordinates[display_index[mask], 1], s=1.2,
                           color=LINEAGE_COLORS[group], alpha=0.56 if group != "Endothelial" else 0.8,
                           linewidths=0, label=group, rasterized=True)
    axis_a.legend(frameon=False, markerscale=5, fontsize=8, ncol=2, loc="upper right")
    axis_a.set(title="Independent PDAC cohort: author-provided major lineages", xlabel="UMAP 1", ylabel="UMAP 2")
    axis_a.set_xticks([]); axis_a.set_yticks([])
    panel_label(axis_a, "A")

    axis_b = figure.add_subplot(grid[0, 1])
    background = np.random.choice(np.arange(data.n_obs), size=min(90000, data.n_obs), replace=False)
    axis_b.scatter(coordinates[background, 0], coordinates[background, 1], s=1.0, color="#D9D9D9", alpha=0.22, linewidths=0, rasterized=True)
    positive = epas > 0
    order = np.flatnonzero(positive)
    values = np.clip(epas[order], 0, np.quantile(epas[positive], 0.995))
    cmap = LinearSegmentedColormap.from_list("epas", ["#FEE8C8", "#E34A33", "#8B0000"])
    dots = axis_b.scatter(coordinates[order, 0], coordinates[order, 1], c=values, cmap=cmap, s=1.8, alpha=0.68, linewidths=0, rasterized=True)
    colorbar = figure.colorbar(dots, ax=axis_b, fraction=0.045, pad=0.02)
    colorbar.set_label("EPAS1 normalized expression", fontsize=9)
    axis_b.set(title="EPAS1 feature map", xlabel="UMAP 1", ylabel="UMAP 2")
    axis_b.set_xticks([]); axis_b.set_yticks([])
    panel_label(axis_b, "B")

    axis_c = figure.add_subplot(grid[1, 0])
    paired = patient_lineage.pivot(index="pid", columns="lineage", values="EPAS1_mean").dropna(subset=["Endothelial"])
    other = paired.drop(columns=["Endothelial"], errors="ignore").median(axis=1)
    paired = pd.DataFrame({"other": other, "endothelial": paired["Endothelial"]}).dropna().sort_values("endothelial")
    jitter = np.linspace(-0.04, 0.04, len(paired))
    for position, (_, row) in enumerate(paired.iterrows()):
        axis_c.plot([0 + jitter[position], 1 + jitter[position]], [row["other"], row["endothelial"]], color="#BEBEBE", linewidth=0.65, zorder=1)
    axis_c.scatter(np.zeros(len(paired)) + jitter, paired["other"], color="#9E9E9E", s=24, zorder=2, label="Patient median of other lineages")
    axis_c.scatter(np.ones(len(paired)) + jitter, paired["endothelial"], color="#177E89", s=28, zorder=3, label="Endothelial")
    axis_c.set_xticks([0, 1], ["Other lineages", "Endothelial"])
    axis_c.set_ylabel("Patient-level mean EPAS1 expression")
    axis_c.set_title(f"Paired lineage localization (n={int(summary.paired_n_patients)}; P={summary.paired_wilcoxon_greater_p:.2e})")
    axis_c.legend(frameon=False, fontsize=8, loc="upper left")
    axis_c.spines[["top", "right"]].set_visible(False)
    panel_label(axis_c, "C")

    axis_d = figure.add_subplot(grid[1, 1])
    subtype = subtype.sort_values("EPAS1_mean", ascending=True)
    size = 50 + 360 * subtype["EPAS1_detection"].to_numpy()
    scatter = axis_d.scatter(subtype["EPAS1_mean"], np.arange(len(subtype)), s=size, c=subtype["ANGPT2_mean"], cmap="YlOrRd", edgecolor="white", linewidth=0.8)
    axis_d.set_yticks(np.arange(len(subtype)), subtype["endothelial_subtype"].str.replace("Endothelial-Vascular_", "", regex=False), fontsize=9)
    axis_d.set_xlabel("Mean EPAS1 expression")
    axis_d.set_title("Source-study endothelial-state labels")
    axis_d.grid(axis="x", color="#E5E5E5", linewidth=0.7)
    axis_d.spines[["top", "right", "left"]].set_visible(False)
    colorbar = figure.colorbar(scatter, ax=axis_d, fraction=0.045, pad=0.02)
    colorbar.set_label("Mean ANGPT2 expression", fontsize=9)
    axis_d.text(0.02, -0.18, "Point area: EPAS1 detection fraction", transform=axis_d.transAxes, fontsize=8, color="#4D4D4D")
    panel_label(axis_d, "D")

    figure.suptitle("Independent PDAC cohort supports endothelial enrichment of EPAS1", fontsize=16, fontweight="bold")
    for extension in ["pdf", "png"]:
        figure.savefig(OUT / f"External_Figure_GSE202051_endothelial_localization.{extension}", dpi=450, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
