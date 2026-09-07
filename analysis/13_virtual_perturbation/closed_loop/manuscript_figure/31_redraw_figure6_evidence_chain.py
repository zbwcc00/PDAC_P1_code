from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BASE = Path(r"D:/PDAC_P1")
OUT = BASE / "analysis/13_virtual_perturbation/closed_loop/manuscript_figure"


def card(ax, x, y, width, height, tag, heading, body, colour):
    box = FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.018",
        transform=ax.transAxes, facecolor="#F8FAFC", edgecolor=colour, linewidth=2.0,
    )
    ax.add_patch(box)
    ax.text(x + 0.018, y + height - 0.035, tag, transform=ax.transAxes,
            fontsize=12, fontweight="bold", color=colour, va="top")
    ax.text(x + 0.065, y + height - 0.036, heading, transform=ax.transAxes,
            fontsize=11, fontweight="bold", color="#172B4D", va="top")
    ax.text(x + width / 2, y + height / 2 - 0.005, body, transform=ax.transAxes,
            fontsize=9.5, color="#243B53", ha="center", va="center", linespacing=1.35)


def arrow(ax, start, end, colour="#607D8B", rad=0.0):
    ax.add_patch(FancyArrowPatch(start, end, transform=ax.transAxes,
                                 arrowstyle="-|>", mutation_scale=17, linewidth=1.8,
                                 color=colour, connectionstyle=f"arc3,rad={rad}"))


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "savefig.dpi": 450})
    fig, ax = plt.subplots(figsize=(14, 8.4))
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.10, top=0.87)

    fig.text(0.5, 0.955, "Figure 6. EPAS1-centered evidence chain and therapeutic hypothesis in PDAC",
             ha="center", va="top", fontsize=19, fontweight="bold", color="#102A43")
    fig.text(0.5, 0.915, "A cross-cohort endothelial program is spatially anchored, computationally perturbable, and used to prioritize compounds for testing",
             ha="center", va="top", fontsize=10.5, color="#52606D")

    card(ax, 0.04, 0.58, 0.25, 0.22, "A", "Endothelial discovery",
         "Three independent scRNA-seq cohorts\nconverge on an EPAS1-centered\nvascular program", "#1F78A8")
    card(ax, 0.375, 0.58, 0.25, 0.22, "B", "Spatial anchoring",
         "EPAS1 co-localizes with vascular\nmarkers in image-resolved PDAC\n(GSE282302: 14/14; GSE297144: 8/8)", "#00897B")
    card(ax, 0.71, 0.58, 0.25, 0.22, "C", "Secreted / immune context",
         "Five SecAct candidates replicate;\nANGPT2 communication and\nmyeloid/T-NK neighborhoods are positive", "#2A9D8F")
    card(ax, 0.14, 0.22, 0.25, 0.22, "D", "Virtual perturbation",
         "EPAS1 KO changes endothelial and\nangiogenesis programs in 3/3 cohorts;\nOE remains supportive (1/3)", "#C23B4A")
    card(ax, 0.45, 0.22, 0.25, 0.22, "E", "DrugReflector prioritization",
         "Y-39983: primary candidate\ntriclabendazole: sensitivity candidate\n(recurrent reverse signatures)", "#D97706")
    card(ax, 0.76, 0.22, 0.20, 0.22, "F", "Structure support",
         "Five EPAS1 PAS-B conformers\n5/5 redocking QC pass\nstructure-guided hypothesis", "#6B4C9A")

    arrow(ax, (0.29, 0.69), (0.375, 0.69))
    arrow(ax, (0.625, 0.69), (0.71, 0.69))
    arrow(ax, (0.83, 0.58), (0.30, 0.44), rad=-0.16)
    arrow(ax, (0.39, 0.33), (0.45, 0.33))
    arrow(ax, (0.70, 0.33), (0.76, 0.33))

    ax.text(0.5, 0.505, "EPAS1-centered vascular–immune remodeling hypothesis",
            transform=ax.transAxes, ha="center", va="center", fontsize=11,
            fontweight="bold", color="#16324F")
    ax.plot([0.34, 0.66], [0.505, 0.505], transform=ax.transAxes, color="#9FB3C8", linewidth=1.0)

    ax.text(0.5, 0.075,
            "Evidence grading: observational localization → model-based signaling → computational perturbation → compound prioritization → structural hypothesis",
            transform=ax.transAxes, ha="center", va="center", fontsize=9.2, color="#52606D")
    ax.text(0.5, 0.035,
            "Scope boundary: these analyses do not establish direct causality, target engagement, pharmacologic efficacy, or clinical benefit; experimental validation is required.",
            transform=ax.transAxes, ha="center", va="center", fontsize=9.2, color="#B42318", fontweight="bold")

    out = OUT / "Figure_6_EPAS1_evidence_chain"
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", pad_inches=0.12)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.12)
    print(out.with_suffix(".png"))
    print(out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
