from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Ellipse, Rectangle, Arc

BASE = Path(r"D:/PDAC_P1")
OUT = BASE / "analysis/13_virtual_perturbation/closed_loop/manuscript_figure"
NAVY, TEXT, MUTED = "#17324D", "#29465F", "#64788A"
PALE, GRID = "#F5F8FA", "#DCE6EC"
TEAL, BLUE, CORAL, AMBER, PURPLE = "#168A86", "#2379A8", "#C83D50", "#D98200", "#6B4C9A"

def rounded(ax, x, y, w, h, edge, face=PALE, lw=2.0, radius=0.018):
    ax.add_patch(FancyBboxPatch((x, y), w, h, transform=ax.transAxes,
        boxstyle=f"round,pad=0.008,rounding_size={radius}", facecolor=face,
        edgecolor=edge, linewidth=lw))

def arrow(ax, start, end, colour="#78909C", rad=0.0, lw=1.6, scale=13):
    ax.add_patch(FancyArrowPatch(start, end, transform=ax.transAxes,
        arrowstyle="-|>", mutation_scale=scale, linewidth=lw, color=colour,
        connectionstyle=f"arc3,rad={rad}"))

def txt(ax, x, y, text, size=8.5, colour=TEXT, weight="normal", ha="left", va="center"):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, color=colour,
            fontweight=weight, ha=ha, va=va)

def cohort_glyph(ax, x, y, colour):
    for i, yy in enumerate([y + 0.044, y + 0.020, y - 0.004]):
        ax.add_patch(Circle((x + 0.020, yy), 0.009, transform=ax.transAxes,
            facecolor=colour, edgecolor="white", linewidth=0.8))
        ax.add_patch(Ellipse((x + 0.058, yy), 0.048, 0.016, transform=ax.transAxes,
            facecolor="#DCECF2", edgecolor=colour, linewidth=0.9))
        txt(ax, x + 0.092, yy, f"cohort {i + 1}", size=6.4, colour=MUTED)
    ax.add_patch(Arc((x + 0.185, y + 0.060), 0.040, 0.040, theta1=210, theta2=330,
        transform=ax.transAxes, color=colour, linewidth=1.1))

def spatial_glyph(ax, x, y, colour):
    for r in range(4):
        for c in range(6):
            fc = ["#EAF2F4", "#CBE4E4", "#8CCBC2", colour][(r + c) % 4]
            ax.add_patch(Rectangle((x + c * 0.018, y + r * 0.018), 0.015, 0.015,
                transform=ax.transAxes, facecolor=fc, edgecolor="white", linewidth=0.3))
    ax.add_patch(Circle((x + 0.046, y + 0.045), 0.018, transform=ax.transAxes,
        facecolor="#F4A261", edgecolor="white", linewidth=0.8))
    txt(ax, x, y - 0.019, "EPAS1 + vascular markers", size=6.2, colour=MUTED)

def immune_glyph(ax, x, y, _colour=None):
    ax.add_patch(Ellipse((x + 0.035, y + 0.034), 0.060, 0.038, transform=ax.transAxes,
        facecolor="#D7F0ED", edgecolor=TEAL, linewidth=1.0))
    ax.add_patch(Circle((x + 0.100, y + 0.046), 0.019, transform=ax.transAxes,
        facecolor="#F6D8C8", edgecolor=CORAL, linewidth=1.0))
    ax.add_patch(Circle((x + 0.137, y + 0.018), 0.014, transform=ax.transAxes,
        facecolor="#E9DFF2", edgecolor=PURPLE, linewidth=1.0))
    arrow(ax, (x + 0.065, y + 0.035), (x + 0.088, y + 0.045), colour=AMBER, lw=1.0, scale=8)
    arrow(ax, (x + 0.065, y + 0.027), (x + 0.122, y + 0.021), colour=CORAL, lw=1.0, scale=8)
    txt(ax, x, y - 0.010, "ANGPT2 / myeloid / T–NK", size=6.2, colour=MUTED)

def ko_glyph(ax, x, y, _colour=None):
    for i, h in enumerate([0.030, 0.060, 0.086]):
        ax.add_patch(Rectangle((x + i * 0.035, y), 0.020, h, transform=ax.transAxes,
            facecolor=["#B9CBD7", CORAL, CORAL][i], edgecolor="none"))
    txt(ax, x, y - 0.016, "KO response", size=6.2, colour=MUTED)
    txt(ax, x + 0.095, y + 0.042, "3/3", size=9, colour=CORAL, weight="bold")

def drug_glyph(ax, x, y, _colour=None):
    for yy, col, name in [(y + 0.043, AMBER, "Y-39983"), (y + 0.006, "#E4A95B", "triclabendazole")]:
        ax.add_patch(FancyBboxPatch((x, yy), 0.095, 0.023, transform=ax.transAxes,
            boxstyle="round,pad=0.004,rounding_size=0.008", facecolor=col,
            edgecolor="white", linewidth=0.7))
        txt(ax, x + 0.047, yy + 0.011, name, size=5.8, colour="white", weight="bold", ha="center")
        arrow(ax, (x + 0.104, yy + 0.011), (x + 0.146, yy + 0.011), colour=col, lw=1.0, scale=8)
    txt(ax, x, y - 0.013, "recurrent reverse signatures", size=6.2, colour=MUTED)

def structure_glyph(ax, x, y, _colour=None):
    ax.add_patch(Ellipse((x + 0.050, y + 0.034), 0.105, 0.058, transform=ax.transAxes,
        facecolor="#E8E0F2", edgecolor=PURPLE, linewidth=1.2))
    ax.add_patch(Ellipse((x + 0.075, y + 0.036), 0.045, 0.024, transform=ax.transAxes,
        facecolor="#FFFFFF", edgecolor=PURPLE, linewidth=1.0))
    ax.add_patch(Circle((x + 0.075, y + 0.036), 0.007, transform=ax.transAxes,
        facecolor=AMBER, edgecolor="white", linewidth=0.7))
    txt(ax, x, y - 0.012, "PAS-B pocket · 5/5 QC", size=6.2, colour=MUTED)

def card(ax, x, y, w, h, tag, heading, body, colour, glyph):
    rounded(ax, x, y, w, h, colour)
    txt(ax, x + 0.018, y + h - 0.030, tag, size=11, colour=colour, weight="bold", va="top")
    txt(ax, x + 0.060, y + h - 0.030, heading, size=10.2, colour=NAVY, weight="bold", va="top")
    glyph(ax, x + 0.020, y + 0.050, colour)
    txt(ax, x + 0.175, y + h / 2 + 0.010, body, size=6.8, colour=TEXT, va="center")

def central_niche(ax):
    rounded(ax, 0.365, 0.405, 0.270, 0.140, "#9AAEBC", face="#FBFCFD", lw=1.4, radius=0.022)
    txt(ax, 0.500, 0.523, "EPAS1-centered vascular niche", size=9.5, colour=NAVY, weight="bold", ha="center")
    ax.add_patch(Ellipse((0.455, 0.462), 0.092, 0.044, transform=ax.transAxes,
        facecolor="#D6EEF0", edgecolor=TEAL, linewidth=1.6))
    txt(ax, 0.455, 0.462, "EPAS1", size=8, colour=TEAL, weight="bold", ha="center")
    for xx, yy, fc, ec in [(0.410, 0.437, "#E8DDF0", PURPLE), (0.565, 0.454, "#F6D8C8", CORAL), (0.575, 0.490, "#F9E6B9", AMBER)]:
        ax.add_patch(Circle((xx, yy), 0.015, transform=ax.transAxes, facecolor=fc, edgecolor=ec, linewidth=1.0))
    arrow(ax, (0.490, 0.462), (0.550, 0.454), colour=CORAL, lw=1.1, scale=9)
    arrow(ax, (0.485, 0.469), (0.565, 0.490), colour=AMBER, lw=1.1, scale=9)
    arrow(ax, (0.443, 0.453), (0.420, 0.440), colour=PURPLE, lw=1.0, scale=8)
    txt(ax, 0.500, 0.421, "vascular state ↔ secreted / immune context", size=6.6, colour=MUTED, ha="center")

def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "savefig.dpi": 600,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(14.5, 8.2))
    ax.axis("off")
    fig.subplots_adjust(left=0.015, right=0.985, bottom=0.055, top=0.86)
    fig.text(0.5, 0.955, "Figure 6. EPAS1-centered evidence chain and therapeutic hypothesis in PDAC",
        ha="center", va="top", fontsize=18, fontweight="bold", color="#102A43")
    fig.text(0.5, 0.918, "Replicated endothelial state → spatial niche → perturbation → compound prioritization → structure-guided hypothesis",
        ha="center", va="top", fontsize=10.2, color="#52606D")
    card(ax, 0.035, 0.605, 0.285, 0.235, "A", "Endothelial discovery", "3 scRNA-seq cohorts\nconverge on EPAS1\nvascular program", BLUE, cohort_glyph)
    card(ax, 0.357, 0.605, 0.285, 0.235, "B", "Spatial anchoring", "EPAS1 co-localizes\nwith vascular markers\n14/14 + 8/8 patients", TEAL, spatial_glyph)
    card(ax, 0.679, 0.605, 0.285, 0.235, "C", "Secreted / immune context", "5 SecAct candidates\nANGPT2 + myeloid/T–NK\nneighborhoods positive", TEAL, immune_glyph)
    card(ax, 0.100, 0.105, 0.250, 0.235, "D", "Virtual perturbation", "KO shifts endothelial\nangiogenesis programs\n3/3 cohorts; OE 1/3", CORAL, ko_glyph)
    card(ax, 0.370, 0.105, 0.250, 0.235, "E", "DrugReflector prioritization", "Y-39983 primary;\ntriclabendazole\nsensitivity candidate", AMBER, drug_glyph)
    card(ax, 0.679, 0.105, 0.285, 0.235, "F", "Structure support", "5 PAS-B conformers;\n5/5 redocking pass;\nstructure hypothesis", PURPLE, structure_glyph)
    central_niche(ax)
    arrow(ax, (0.320, 0.855), (0.357, 0.855), colour=BLUE, lw=1.7)
    arrow(ax, (0.642, 0.855), (0.679, 0.855), colour=TEAL, lw=1.7)
    arrow(ax, (0.825, 0.605), (0.635, 0.510), colour=TEAL, rad=0.12, lw=1.4)
    arrow(ax, (0.365, 0.460), (0.285, 0.340), colour=CORAL, rad=0.08, lw=1.4)
    arrow(ax, (0.350, 0.345), (0.370, 0.345), colour=CORAL, lw=1.7)
    arrow(ax, (0.620, 0.345), (0.679, 0.345), colour=AMBER, lw=1.7)
    arrow(ax, (0.500, 0.405), (0.500, 0.340), colour="#91A4B2", lw=1.0, scale=9)
    txt(ax, 0.500, 0.580, "evidence integration", size=8.4, colour=NAVY, weight="bold", ha="center")
    ax.plot([0.415, 0.585], [0.570, 0.570], transform=ax.transAxes, color=GRID, linewidth=1.0)
    txt(ax, 0.500, 0.072, "Evidence grading: observational localization → model-based signaling → computational perturbation → compound prioritization → structural hypothesis", size=8.5, colour=MUTED, ha="center")
    txt(ax, 0.500, 0.033, "Scope boundary: these analyses do not establish direct causality, target engagement, pharmacologic efficacy, or clinical benefit; experimental validation is required.", size=8.4, colour="#B42318", weight="bold", ha="center")
    out = OUT / "Figure_6_EPAS1_evidence_chain"
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", pad_inches=0.12)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.12)
    print(out.with_suffix(".png"))
    print(out.with_suffix(".pdf"))

if __name__ == "__main__":
    main()
