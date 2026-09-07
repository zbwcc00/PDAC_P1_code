from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image

BASE = Path(r"D:/PDAC_P1")
SPATIAL = BASE / "analysis/09_spatial_validation"
CLOSED = BASE / "analysis/13_virtual_perturbation/closed_loop"
DOCK = BASE / "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16"
OUT = CLOSED / "manuscript_figure"

def heatmap(ax, data, cmap, vmin, vmax, labels, title, cbar_label=None):
    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(data.shape[1]), labels[1])
    ax.set_yticks(range(data.shape[0]), labels[0])
    ax.tick_params(axis="both", labelsize=7)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if np.isfinite(data[i, j]):
                ax.text(j, i, f"{data[i, j]:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    if cbar_label:
        cbar = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
        cbar.ax.tick_params(labelsize=7)
        cbar.set_label(cbar_label, fontsize=7)

def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 450})
    fig = plt.figure(figsize=(15, 13), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 0.95], width_ratios=[1, 1])

    # A: spatial anchor image, retained as the empirical localization panel.
    ax = fig.add_subplot(grid[0, 0])
    anchor = Image.open(SPATIAL / "figures/Figure_main_spatial_EPAS1_vascular_anchor.png")
    ax.imshow(anchor)
    ax.axis("off")
    ax.set_title("A  Independent spatial EPAS1–vascular anchoring", loc="left", fontsize=10, fontweight="bold")

    # B: SecAct cross-cohort patient-level replication.
    ax = fig.add_subplot(grid[0, 1])
    sec = pd.read_csv(SPATIAL / "SecAct/SecAct_five_targets_patient_summary_1000.tsv", sep="\t")
    core = ["SPARCL1", "VWF", "MMRN2", "ANGPT2", "IL33"]
    sec = sec[sec["secreted"].isin(core)]
    tab = sec.pivot(index="secreted", columns="dataset", values="positive_fraction").reindex(core)
    heatmap(ax, tab.values, "Blues", 0, 1, (tab.index, tab.columns),
            "B  Patient-level SecAct replication (1000 permutations)", "Positive fraction")
    ax.set_xlabel("Both cohorts show 8/8 or 14/14 positive patients", fontsize=7)

    # C: virtual EPAS1 perturbation in endothelial cells.
    ax = fig.add_subplot(grid[1, 0])
    mech = pd.read_csv(CLOSED / "epas1_mechanism_program_summary.tsv", sep="\t")
    mech = mech[(mech["cell_type"] == "endothelial") & (mech["program"].isin(["endothelial", "angiogenesis", "TNF_NFkB", "myeloid_inflammation", "CAF_ECM", "cytotoxic_TNK", "T_exhaustion"]))]
    programs = ["endothelial", "angiogenesis", "TNF_NFkB", "myeloid_inflammation", "CAF_ECM", "cytotoxic_TNK", "T_exhaustion"]
    z = np.array([[float(mech[(mech.program == p) & (mech.perturbation == q)].iloc[0].median_abs_Z) for q in ["KO", "OE"]] for p in programs])
    heatmap(ax, z, "YlOrRd", 0, max(1.8, float(np.nanmax(z))), (programs, ["KO", "OE"]),
            "C  Endothelial EPAS1 virtual perturbation", "Median |Z|")
    ax.set_xlabel("Multi-dataset computational response; strongest reproducibility in KO endothelial/angiogenesis programs", fontsize=7)

    # D: immune ecology context.
    ax = fig.add_subplot(grid[1, 1])
    eco = pd.read_csv(SPATIAL / "SecAct/SecAct_spatial_ecology_context_summary.tsv", sep="\t").iloc[0]
    vals = [float(eco["median_myeloid_delta"]), float(eco["median_T_NK_delta"])]
    bars = ax.bar(["Myeloid", "T/NK"], vals, color=["#0072B2", "#009E73"], width=0.55)
    ax.axhline(0, color="#555", lw=0.8)
    ax.set_ylabel("EPAS1-high minus EPAS1-low score")
    ax.set_title("D  Immune neighborhood context in GSE282302", loc="left", fontsize=10, fontweight="bold")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.003, f"{val:.3f}", ha="center", fontsize=9)
    ax.text(0.03, 0.92, "14/14 patients positive for both programs", transform=ax.transAxes, fontsize=8)

    # E: structure-supported compound prioritization.
    ax = fig.add_subplot(grid[2, 0])
    dock = pd.read_csv(DOCK / "final_exhaustiveness16_summary.tsv", sep="\t")
    y = [dock.loc[dock["compound"] == cid, "affinity_best_kcal_mol"].astype(float).tolist() for cid in ["BRD-K56751279", "BRD-K81916719"]]
    ax.boxplot(y, positions=[1, 2], widths=0.45, patch_artist=True,
               boxprops={"facecolor": "#E6E6E6"}, medianprops={"color": "#D55E00", "linewidth": 2})
    for pos, values in zip([1, 2], y):
        ax.scatter(np.full(len(values), pos), values, color="#333", s=22, zorder=3)
    ax.set_xticks([1, 2], ["Y-39983", "Triclabendazole"])
    ax.set_ylabel("Best-mode Vina affinity (kcal/mol)")
    ax.set_title("E  DrugReflector-to-EPAS1 PAS-B docking", loc="left", fontsize=10, fontweight="bold")
    ax.text(0.02, 0.05, "5 conformers; exhaustiveness=16; 5/5 co-crystal RMSD QC passed", transform=ax.transAxes, fontsize=7)

    # F: explicit evidence chain and scope boundary.
    ax = fig.add_subplot(grid[2, 1])
    ax.axis("off")
    ax.set_title("F  EPAS1 closed-loop interpretation", loc="left", fontsize=10, fontweight="bold")
    boxes = [(0.03, 0.58, 0.26, 0.24, "Spatial anchor\nEPAS1 near vascular markers"),
             (0.37, 0.58, 0.26, 0.24, "SecAct + ecology\nsecreted signals + immune context"),
             (0.71, 0.58, 0.26, 0.24, "Virtual perturbation\nendothelial/angiogenesis response"),
             (0.24, 0.12, 0.26, 0.24, "DrugReflector\nY-39983 primary; triclabendazole sensitivity"),
             (0.58, 0.12, 0.26, 0.24, "Docking\nstructure-supported hypothesis")]
    for x, y0, w, h, text in boxes:
        patch = FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.015", facecolor="#F4F7FA", edgecolor="#4C78A8", linewidth=1.2, transform=ax.transAxes)
        ax.add_patch(patch)
        ax.text(x + w/2, y0 + h/2, text, ha="center", va="center", fontsize=8, transform=ax.transAxes)
    arrows = [((0.29, 0.70), (0.37, 0.70)), ((0.63, 0.70), (0.71, 0.70)), ((0.84, 0.58), (0.69, 0.36)), ((0.50, 0.58), (0.39, 0.36))]
    for start, end in arrows:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12, linewidth=1.1, color="#777", transform=ax.transAxes))
    ax.text(0.5, 0.02, "Computational prioritization only; no direct binding, efficacy, or causality claim", ha="center", fontsize=7.5, color="#8B0000", transform=ax.transAxes)

    fig.suptitle("EPAS1-centered vascular–immune signaling and therapeutic prioritization in PDAC", fontsize=15, fontweight="bold")
    out = OUT / "EPAS1_closed_loop_main_mechanism_figure"
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", pad_inches=0.1)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.1)
    print(out.with_suffix('.png'))
    print(out.with_suffix('.pdf'))

if __name__ == "__main__":
    main()
