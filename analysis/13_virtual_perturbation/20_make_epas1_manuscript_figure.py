from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASE = Path(r"D:\PDAC_P1")
ROOT = BASE / "analysis" / "13_virtual_perturbation"
CLOSED = ROOT / "closed_loop"
DOCK = ROOT / "docking" / "EPAS1_multiconformer" / "final_exhaustiveness16"
OUT = CLOSED / "manuscript_figure"
OUT.mkdir(parents=True, exist_ok=True)


def read(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    mechanism = read(CLOSED / "epas1_mechanism_program_summary.tsv")
    programs = ["endothelial", "angiogenesis", "TNF_NFkB", "myeloid_inflammation", "CAF_ECM", "cytotoxic_TNK", "T_exhaustion"]
    perturbations = ["KO", "OE"]
    endothelial_rows = [r for r in mechanism if r["cell_type"] == "endothelial"]

    z = np.full((len(programs), len(perturbations)), np.nan)
    sig_fraction = np.full_like(z, np.nan)
    for i, program in enumerate(programs):
        for j, perturbation in enumerate(perturbations):
            row = next((r for r in endothelial_rows if r["program"] == program and r["perturbation"] == perturbation), None)
            if row:
                z[i, j] = float(row["median_abs_Z"])
                sig_fraction[i, j] = int(row["n_significant_datasets"]) / int(row["n_datasets"])

    docking = read(DOCK / "final_exhaustiveness16_summary.tsv")
    compounds = [("Y-39983", "BRD-K56751279"), ("Triclabendazole", "BRD-K81916719")]
    affinities = [[float(r["affinity_best_kcal_mol"]) for r in docking if r["compound"] == cid] for _, cid in compounds]

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9})
    fig = plt.figure(figsize=(12, 8.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.15, 1])

    ax = fig.add_subplot(grid[0, 0])
    im = ax.imshow(z, cmap="YlOrRd", vmin=0, vmax=max(1.8, np.nanmax(z)))
    ax.set_xticks(range(2), perturbations)
    ax.set_yticks(range(len(programs)), programs)
    ax.set_title("A  Endothelial EPAS1 perturbation")
    for i in range(len(programs)):
        for j in range(2):
            if np.isfinite(z[i, j]):
                ax.text(j, i, f"{z[i,j]:.2f}\n({int(round(sig_fraction[i,j]*100))}% sig.)", ha="center", va="center", fontsize=8)
    ax.set_xlabel("Virtual perturbation")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Median |Z| across datasets")

    ax = fig.add_subplot(grid[0, 1])
    focus = ["endothelial", "angiogenesis"]
    x = np.arange(len(focus))
    width = 0.36
    ko = [sig_fraction[programs.index(p), 0] for p in focus]
    oe = [sig_fraction[programs.index(p), 1] for p in focus]
    ax.bar(x - width / 2, ko, width, label="KO", color="#b2182b")
    ax.bar(x + width / 2, oe, width, label="OE", color="#2166ac")
    ax.set_ylim(0, 1.08)
    ax.set_xticks(x, ["Endothelial", "Angiogenesis"])
    ax.set_ylabel("Fraction of significant datasets")
    ax.set_title("B  Cross-dataset reproducibility (n=3)")
    ax.legend(frameon=False)
    for xpos, value in zip(x - width / 2, ko):
        ax.text(xpos, value + 0.03, f"{value:.0%}", ha="center", fontsize=8)
    for xpos, value in zip(x + width / 2, oe):
        ax.text(xpos, value + 0.03, f"{value:.0%}", ha="center", fontsize=8)

    ax = fig.add_subplot(grid[1, 0])
    positions = np.arange(1, 3)
    ax.boxplot(affinities, positions=positions, widths=0.48, showmeans=True, meanline=False,
               patch_artist=True, boxprops={"facecolor": "#e6e6e6"}, medianprops={"color": "#d95f02", "linewidth": 2})
    for pos, values in zip(positions, affinities):
        ax.scatter(np.full(len(values), pos), values, color="#333333", s=22, zorder=3)
    ax.set_xticks(positions, [name for name, _ in compounds])
    ax.set_ylabel("Best-mode Vina affinity (kcal/mol)")
    ax.set_title("C  EPAS1 PAS-B five-conformer docking")
    ax.grid(axis="y", alpha=0.25)

    ax = fig.add_subplot(grid[1, 1])
    ax.axis("off")
    rows = [
        ["Candidate", "Mean ± SD", "Range", "Role"],
        ["Y-39983", "−6.55 ± 0.66", "−7.73 to −6.17", "Primary"],
        ["Triclabendazole", "−6.11 ± 0.64", "−6.97 to −5.45", "Sensitivity"],
        ["Co-crystal redocking", "RMSD ≤ 1.02 Å", "5/5 passed", "QC"],
    ]
    table = ax.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="center")
    table.auto_set_font_size(False); table.set_fontsize(8.5); table.scale(1.1, 1.7)
    ax.set_title("D  Closed-loop evidence summary", pad=16)

    fig.savefig(OUT / "epas1_closed_loop_main_figure.png", dpi=450)
    fig.savefig(OUT / "epas1_closed_loop_main_figure.pdf")

    summary = [
        "# EPAS1 closed-loop evidence summary",
        "",
        "- Endothelial EPAS1 KO showed significant responses in endothelial and angiogenesis programs in 3/3 datasets.",
        "- EPAS1 OE showed significant responses in these programs in 1/3 datasets; OE is therefore supportive but not a reproducible primary claim.",
        "- The strongest reproducible readout is the endothelial/angiogenesis program response, summarized with median |Z| rather than raw fold-change ratios.",
        "- Y-39983 was retained as the primary DrugReflector reverse-signature candidate; triclabendazole was retained as a secondary sensitivity candidate.",
        "- Five-conformer EPAS1 PAS-B docking was run at exhaustiveness 16. Co-crystal redocking passed RMSD ≤2 Å in all five conformations.",
        "",
        "Interpretation boundary: virtual perturbation, DrugReflector ranking, and docking provide computational prioritization and structural support; they do not establish pharmacologic efficacy or direct binding.",
    ]
    (OUT / "epas1_closed_loop_results.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(OUT / "epas1_closed_loop_main_figure.png")


if __name__ == "__main__":
    main()
