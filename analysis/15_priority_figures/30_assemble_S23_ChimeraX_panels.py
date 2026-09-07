from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image


out_dir = Path(r"D:/PDAC_P1/analysis/15_priority_figures")
y39983 = Image.open(out_dir / "EPAS1_Y39983_ChimeraX.png")
triclabendazole = Image.open(out_dir / "EPAS1_triclabendazole_ChimeraX.png")

fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.25), dpi=300)
for axis, image in zip(axes, [y39983, triclabendazole]):
    axis.imshow(image)
    axis.axis("off")
axes[0].set_title("Y-39983 | EPAS1 PAS-B (6D0B) | Vina −7.729 kcal/mol", fontsize=11, fontweight="bold", pad=8)
axes[1].set_title("Triclabendazole | EPAS1 PAS-B (6CZW) | Vina −6.967 kcal/mol", fontsize=11, fontweight="bold", pad=8)
fig.suptitle("Supplementary Figure S23. ChimeraX representative EPAS1–ligand complexes", fontsize=15, fontweight="bold", y=0.99)
fig.text(0.5, 0.025, "Protein: navy ribbon; ligand: cyan (Y-39983) or magenta (triclabendazole); tan sticks: nearby pocket residues. Dashed cyan lines are ChimeraX geometry-based hydrogen-bond candidates after protonation; absence is retained as a result.", ha="center", va="bottom", fontsize=8.3, color="#444444")
fig.subplots_adjust(left=0.01, right=0.99, top=0.87, bottom=0.09, wspace=0.02)
fig.savefig(out_dir / "Figure_S23_EPAS1_complex_binding_modes.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_dir / "Figure_S23_EPAS1_complex_binding_modes.pdf", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
