from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(r"D:/第二篇大论文")
OUT = BASE / "analysis/15_priority_figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.labelsize": 9, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})
BLUE = "#0072B2"
ORANGE = "#D55E00"
TEAL = "#009E73"
GREY = "#8C8C8C"


def save(fig: plt.Figure, stem: str):
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    plt.close(fig)


def discovery_overview():
    x = pd.read_csv(BASE / "analysis/02_scrna_qc/PDAC_scRNA_dataset_summary.tsv", sep="\t")
    comp = pd.read_csv(BASE / "analysis/03_scrna_annotation/PDAC_scRNA_patient_lineage_summary.tsv", sep="\t")
    lineages = ["endothelial", "fibroblast_CAF", "myeloid", "T_NK", "B_cell", "epithelial", "malignant_epithelial", "acinar", "mast", "unknown_ambiguous"]
    mean_comp = comp.groupby("dataset")[[f"prop_{z}" for z in lineages]].mean().reindex(x.dataset).fillna(0)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw={"width_ratios": [1.0, 1.35]})
    axes[0].bar(x.dataset, x.samples, color=BLUE, label="Samples")
    twin = axes[0].twinx()
    twin.spines["right"].set_position(("outward", 38))
    twin.plot(x.dataset, x.total_cells / 1000, color=ORANGE, marker="o", lw=2, label="Cells (×1000)")
    axes[0].set_ylabel("Samples")
    twin.set_ylabel("")
    axes[0].set_title("Discovery cohort scale")
    axes[0].tick_params(axis="x", rotation=35)
    for i, row in x.iterrows():
        axes[0].text(i, row.samples + 0.7, f"n={row.samples}", ha="center", fontsize=8)
    colors = [BLUE, "#56B4E9", TEAL, "#CC79A7", "#E69F00", "#999999", "#264653", "#F0E442", "#D55E00", "#BBBBBB"]
    bottom = np.zeros(len(mean_comp))
    for col, color in zip(mean_comp.columns, colors):
        vals = mean_comp[col].to_numpy()
        axes[1].bar(mean_comp.index, vals, bottom=bottom, color=color, label=col.replace("prop_", ""))
        bottom += vals
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Mean patient cell fraction")
    axes[1].set_title("Cellular composition")
    axes[1].tick_params(axis="x", rotation=35)
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=7)
    fig.suptitle("Figure S6. Multi-cohort single-cell discovery overview", y=1.03, fontweight="bold")
    save(fig, "Figure_S6_scRNA_discovery_overview")


def convergence_heatmap():
    x = pd.read_csv(BASE / "analysis/07_candidate_screen/candidate_genes/endothelial_all_merged.tsv", sep="\t")
    genes = ["EPAS1", "TACC1", "MARCKS", "HERPUD1"]
    cols = ["edge_logFC_gse155698", "deseq_logFC_gse155698", "edge_logFC_gse212966", "deseq_logFC_gse212966"]
    labels = ["GSE155698\nedgeR", "GSE155698\nDESeq2", "GSE212966\nedgeR", "GSE212966\nDESeq2"]
    mat = x.set_index("gene_id").reindex(genes)[cols].astype(float)
    mat.columns = labels
    fig, ax = plt.subplots(figsize=(6.0, 2.7))
    image = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-5, vmax=5, aspect="auto")
    ax.set_xticks(np.arange(mat.shape[1]), mat.columns)
    ax.set_yticks(np.arange(mat.shape[0]), mat.index)
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            value = mat.iloc[row, col]
            ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, shrink=.8, label="logFC")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Figure S7. Cross-cohort and cross-method endothelial candidate convergence")
    save(fig, "Figure_S7_endothelial_candidate_convergence")


def bulk_validation():
    paired = pd.read_csv(BASE / "analysis/08_program_validation/GSE62452/GSE62452_program_scores.tsv", sep="\t")
    contrasts = pd.read_csv(BASE / "analysis/08_program_validation/GSE71729/GSE71729_program_contrasts.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
    for _, grp in paired[paired.group.isin(["tumor", "adjacent"])].groupby("pair_id"):
        if set(grp.group) == {"tumor", "adjacent"}:
            vals = grp.set_index("group").loc[["adjacent", "tumor"], "endothelial"]
            axes[0].plot([0, 1], vals, color="#BBBBBB", lw=.8, zorder=1)
            axes[0].scatter([0, 1], vals, color=[GREY, BLUE], s=16, zorder=2)
    axes[0].set_xticks([0, 1], ["Adjacent", "Tumor"])
    axes[0].set_ylabel("Endothelial program score")
    axes[0].set_title("GSE62452 paired samples (n=45 pairs)")
    axes[0].text(.5, .97, "paired logFC=0.248; P=1.17×10⁻⁵", transform=axes[0].transAxes, ha="center", va="top", fontsize=8)
    plot = contrasts[contrasts.contrast == "primary_vs_normal"].set_index("program").reindex(["endothelial", "T_NK"])
    bars = axes[1].bar(plot.index, plot.logFC, color=[BLUE, TEAL], width=.6)
    axes[1].axhline(0, color="black", lw=.7)
    axes[1].set_ylabel("Primary vs normal logFC")
    axes[1].set_title("GSE71729 primary vs normal")
    axes[1].tick_params(axis="x", rotation=20)
    for bar, p in zip(bars, plot["P.Value"]):
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + .015, f"P={p:.2g}", ha="center", fontsize=8)
    fig.suptitle("Figure S8. External bulk validation of vascular and immune programs", y=1.04, fontweight="bold")
    save(fig, "Figure_S8_bulk_program_validation")


def spatial_ecology():
    x = pd.read_csv(BASE / "analysis/09_spatial_validation/GSE282302/spatial_ecology/GSE282302_spatial_ecology_patient.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(4.8, 3.1))
    for _, row in x.iterrows():
        ax.plot([0, 1], [row.median_myeloid_delta, row.median_T_NK_delta], color="#BBBBBB", lw=.8, zorder=1)
        ax.scatter([0, 1], [row.median_myeloid_delta, row.median_T_NK_delta], color=[BLUE, ORANGE], s=20, zorder=2)
    ax.axhline(0, color="black", lw=.7)
    ax.set_xticks([0, 1], ["Myeloid", "T/NK"])
    ax.set_ylabel("EPAS1-high − EPAS1-low neighborhood score")
    ax.set_title("Figure S9. Patient-level spatial immune ecology (GSE282302)")
    ax.text(.02, .96, "14/14 patients positive for both programs\nmedian Δ=0.070 and 0.089", transform=ax.transAxes, va="top", fontsize=8)
    save(fig, "Figure_S9_spatial_ecology_patient_level")


def secact_panel():
    x = pd.read_csv(BASE / "analysis/09_spatial_validation/SecAct/SecAct_five_targets_patient_summary_1000.tsv", sep="\t")
    genes = ["SPARCL1", "VWF", "MMRN2", "ANGPT2", "IL33", "IGFBP7"]
    wide = x.pivot(index="secreted", columns="dataset", values="positive_fraction").reindex(genes)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), gridspec_kw={"width_ratios": [1.45, 1]})
    xx = np.arange(len(genes))
    width = .36
    axes[0].bar(xx - width/2, wide["GSE282302"], width, color=BLUE, label="GSE282302 (n=14)")
    axes[0].bar(xx + width/2, wide["GSE297144"], width, color=ORANGE, label="GSE297144 (n=8)")
    axes[0].set_xticks(xx, genes, rotation=35, ha="right")
    axes[0].set_ylim(0, 1.12)
    axes[0].set_ylabel("Patient/sample positive fraction")
    axes[0].set_title("SecAct replication")
    axes[0].legend(frameon=False, fontsize=7, loc="upper left")
    for i, gene in enumerate(genes):
        vals = [wide.loc[gene, ds] for ds in ["GSE282302", "GSE297144"]]
        if not all(value > .95 for value in vals):
            for j, val in enumerate(vals):
                axes[0].text(i + (j-.5)*width, val + .03,
                             f"{val:.0%}", ha="center", va="bottom", fontsize=7,
                             color="black", fontweight="bold")
    edges = pd.DataFrame({"receiver": ["Myeloid", "CD8", "B cell"], "edges": [2, 2, 1]})
    axes[1].barh(edges.receiver, edges.edges, color=TEAL)
    axes[1].set_xlabel("Inferred ANGPT2–ITGA5_ITGB1 edges")
    axes[1].set_title("ANGPT2 immune receivers")
    axes[1].invert_yaxis()
    fig.suptitle("Figure S10. Patient-level SecAct and focused ANGPT2 communication", y=1.04, fontweight="bold")
    save(fig, "Figure_S10_secact_angpt2_replication")


if __name__ == "__main__":
    discovery_overview()
    convergence_heatmap()
    bulk_validation()
    spatial_ecology()
    secact_panel()
    print(OUT)
