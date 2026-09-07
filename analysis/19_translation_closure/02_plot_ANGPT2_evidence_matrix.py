from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("D:/PDAC_P1")
OUT = ROOT / "analysis" / "19_translation_closure"


def main() -> None:
    evidence = pd.read_csv(OUT / "ANGPT2_cross_method_evidence_matrix.tsv", sep="\t")
    edges = pd.read_csv(OUT / "ANGPT2_myeloid_CellChat_dataset_audit.tsv", sep="\t")
    correlations = pd.read_csv(OUT / "ANGPT2_patient_source_receiver_correlation.tsv", sep="\t")

    source = evidence.loc[evidence.evidence_layer.eq("Source expression")].copy()
    receiver = evidence.loc[evidence.evidence_layer.eq("Receiver availability")].copy()
    source["positive"] = source.result.str.extract(r"(\d+)/")[0].astype(int)
    source["total"] = source.result.str.extract(r"/(\d+)")[0].astype(int)
    receiver["positive"] = receiver.result.str.extract(r"(\d+)/")[0].astype(int)
    receiver["total"] = receiver.result.str.extract(r"/(\d+)")[0].astype(int)
    secact = evidence.loc[evidence.evidence_layer.eq("SecAct spatial permutation")].copy()
    secact["positive"] = secact.result.str.extract(r"(\d+)/")[0].astype(int)
    secact["total"] = secact.result.str.extract(r"/(\d+)")[0].astype(int)

    fig = plt.figure(figsize=(13.5, 9.25))
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.05, 1.15],
        left=0.08,
        right=0.98,
        bottom=0.11,
        top=0.84,
        wspace=0.27,
        hspace=0.32,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    labels = ["Endothelial\nANGPT2 detected", "Myeloid ITGA5/ITGB1\njointly detected", "CellChat endothelial →\nmyeloid edge", "SecAct ANGPT2\nspatial activity", "Patient-level\nco-variation"]
    values = [3, 3, int(edges.cellchat_edge_present.sum()), 2, 0]
    totals = [3, 3, len(edges), 2, 3]
    colors = ["#177E89", "#177E89", "#177E89", "#177E89", "#BDBDBD"]
    positions = np.arange(len(labels))
    ax_a.barh(positions, totals, color="#E5E5E5", height=0.56)
    ax_a.barh(positions, values, color=colors, height=0.56)
    for y, value, total in zip(positions, values, totals):
        ax_a.text(total + 0.07, y, f"{value}/{total}", va="center", fontsize=10, fontweight="bold")
    ax_a.set_yticks(positions, labels, fontsize=9)
    ax_a.invert_yaxis(); ax_a.set_xlim(0, 3.6); ax_a.set_xticks([0, 1, 2, 3])
    ax_a.set_xlabel("Supporting cohorts / independently run datasets")
    ax_a.set_title("Cross-method evidence coverage", pad=12)
    ax_a.spines[["top", "right", "left"]].set_visible(False)
    ax_a.text(-0.13, 1.025, "A", transform=ax_a.transAxes, fontsize=15, fontweight="bold")

    ax_b = fig.add_subplot(grid[0, 1])
    source_display = source[["dataset_or_cohort", "positive", "total"]].rename(columns={"positive": "source_positive", "total": "source_total"})
    receiver_display = receiver[["dataset_or_cohort", "positive", "total"]].rename(columns={"positive": "receiver_positive", "total": "receiver_total"})
    display = source_display.merge(receiver_display, on="dataset_or_cohort", how="inner")
    display["source_fraction"] = display.source_positive / display.source_total
    display["receiver_fraction"] = display.receiver_positive / display.receiver_total
    positions = np.arange(len(display)); width = 0.34
    ax_b.bar(positions - width / 2, display.source_fraction, width, color="#177E89", label="Endothelial ANGPT2")
    ax_b.bar(positions + width / 2, display.receiver_fraction, width, color="#6BA292", label="Myeloid ITGA5/ITGB1")
    ax_b.set_xticks(positions, display.dataset_or_cohort, rotation=18, ha="right")
    ax_b.set_ylim(0, 1.13); ax_b.set_ylabel("Patient-group detection fraction")
    ax_b.set_title("Source and receiver availability", pad=12)
    ax_b.legend(frameon=False, fontsize=9, loc="lower left")
    ax_b.spines[["top", "right"]].set_visible(False)
    ax_b.text(-0.12, 1.025, "B", transform=ax_b.transAxes, fontsize=15, fontweight="bold")

    ax_c = fig.add_subplot(grid[1, 0])
    values = np.log10(edges.mean_probability.to_numpy())
    ax_c.scatter(values, np.arange(len(edges)), s=95, color="#177E89", edgecolor="white", linewidth=0.8, zorder=3)
    for y, (value, row) in enumerate(zip(values, edges.itertuples())):
        p_text = "P<1e-4" if row.min_p == 0 else f"P={row.min_p:.1e}"
        ax_c.text(value + 0.03, y, p_text, va="center", fontsize=8)
    ax_c.set_yticks(np.arange(len(edges)), edges.dataset)
    ax_c.invert_yaxis(); ax_c.grid(axis="x", color="#E5E5E5", linewidth=0.7)
    ax_c.set_xlabel("log10 mean CellChat communication probability")
    ax_c.set_title("Endothelial ANGPT2 → myeloid ITGA5/ITGB1", pad=10)
    ax_c.spines[["top", "right", "left"]].set_visible(False)
    ax_c.text(-0.13, 1.025, "C", transform=ax_c.transAxes, fontsize=15, fontweight="bold")

    ax_d = fig.add_subplot(grid[1, 1])
    subset = correlations.loc[correlations.receiver.eq("myeloid")].copy().sort_values("dataset")
    ax_d.axvline(0, color="#777777", linewidth=0.9, linestyle="--")
    colors = np.where(subset.BH_q < 0.05, "#177E89", "#BDBDBD")
    ax_d.scatter(subset.spearman_rho, np.arange(len(subset)), s=90, color=colors, edgecolor="white", linewidth=0.8, zorder=3)
    for y, row in enumerate(subset.itertuples()):
        ax_d.text(row.spearman_rho + 0.04, y, f"q={row.BH_q:.2f}", va="center", fontsize=8)
    ax_d.set_yticks(np.arange(len(subset)), subset.dataset); ax_d.invert_yaxis(); ax_d.set_xlim(-1, 1)
    ax_d.set_xlabel("Patient-level source–receiver Spearman rho")
    ax_d.set_title("Sensitivity: no consistent patient-level co-variation", pad=10)
    ax_d.grid(axis="x", color="#E5E5E5", linewidth=0.7); ax_d.spines[["top", "right", "left"]].set_visible(False)
    ax_d.text(-0.12, 1.025, "D", transform=ax_d.transAxes, fontsize=15, fontweight="bold")

    fig.suptitle(
        "ANGPT2–ITGA5/ITGB1 is a multi-method candidate endothelial–myeloid interface",
        fontsize=15,
        fontweight="bold",
        y=0.965,
    )
    fig.text(0.5, 0.028, "CellChat and expression availability support a candidate interface; spatial SecAct supports ANGPT2 activity in EPAS1-high context. These results do not establish direct signaling or causal immune reprogramming.", ha="center", fontsize=8.5, color="#404040")
    for extension in ["pdf", "png"]:
        fig.savefig(OUT / f"Supplementary_Figure_ANGPT2_cross_method_evidence.{extension}", dpi=450, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
