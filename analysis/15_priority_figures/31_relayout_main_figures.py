from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(r"D:\第二篇大论文")
OUT = ROOT / "analysis" / "15_priority_figures" / "main_figures_relayout"
OUT.mkdir(parents=True, exist_ok=True)

FONT_CANDIDATES = [
    Path(r"C:\Windows\Fonts\arialbd.ttf"),
    Path(r"C:\Windows\Fonts\Arial.ttf"),
]
FONT = next((p for p in FONT_CANDIDATES if p.exists()), None)


def font(size):
    return ImageFont.truetype(str(FONT), size) if FONT else ImageFont.load_default()


def load(rel, crop_top=0):
    image = Image.open(ROOT / rel).convert("RGB")
    if crop_top:
        image = image.crop((0, crop_top, image.width, image.height))
    return image


def fit(img, max_w, max_h):
    scale = min(max_w / img.width, max_h / img.height, 1.0)
    return img.resize((round(img.width * scale), round(img.height * scale)), Image.Resampling.LANCZOS)


def compose(name, title, entries, cols=1):
    margin, gap, label_h = 70, 30, 58
    max_w = 1800
    panel_w = (max_w - 2 * margin - (cols - 1) * gap) // cols
    panel_h = 1100
    rows = (len(entries) + cols - 1) // cols
    canvas = Image.new("RGB", (max_w, margin + 80 + rows * (panel_h + label_h) + (rows - 1) * gap + margin), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text((margin, 25), title, fill="#111111", font=font(42))
    for i, (label, rel, crop_top) in enumerate(entries):
        row, col = divmod(i, cols)
        x = margin + col * (panel_w + gap)
        y = margin + 80 + row * (panel_h + label_h + gap)
        draw.text((x, y), label, fill="#1f3b5b", font=font(34))
        panel = fit(load(rel, crop_top), panel_w, panel_h)
        px = x + (panel_w - panel.width) // 2
        py = y + label_h
        canvas.paste(panel, (px, py))
    png = OUT / f"{name}.png"
    pdf = OUT / f"{name}.pdf"
    canvas.save(png, dpi=(300, 300), optimize=True)
    canvas.save(pdf, "PDF", resolution=300.0)
    print(png)
    print(pdf)


def main():
    compose(
        "Figure_1_scRNA_discovery_convergence",
        "Figure 1. Multi-cohort single-cell discovery of the EPAS1 endothelial axis",
        [
            ("A  Annotation and CopyKAT/CNV malignant-cell assessment", "analysis/15_priority_figures/Figure_S15_scRNA_annotation_CopyKAT_CNV.png", 145),
            ("B  Cross-cohort and cross-method endothelial candidate convergence", "analysis/15_priority_figures/Figure_S7_endothelial_candidate_convergence.png", 115),
        ],
        cols=1,
    )
    compose(
        "Figure_2_bulk_spatial_validation",
        "Figure 2. Bulk and spatial validation of EPAS1 vascular localization",
        [
            ("A  External bulk program validation", "analysis/15_priority_figures/Figure_S8_bulk_program_validation.png", 110),
            ("B  Image-resolved EPAS1/vascular-marker spatial localization", "analysis/15_priority_figures/Figure_S21_spatial_image_montage.png", 65),
        ],
        cols=1,
    )
    compose(
        "Figure_3_secact_immune_ecology",
        "Figure 3. Secreted signaling and immune spatial ecology",
        [
            ("A  SecAct replication and ANGPT2 communication", "analysis/15_priority_figures/Figure_S10_secact_angpt2_replication.png", 115),
            ("B  Patient-level myeloid and T/NK neighborhood associations", "analysis/15_priority_figures/Figure_S9_spatial_ecology_patient_level.png", 105),
        ],
        cols=1,
    )
    compose(
        "Figure_4_virtual_perturbation",
        "Figure 4. Virtual perturbation prioritizes an endothelial EPAS1 program",
        [("A  scTenifoldKnk and scTenifoldNet program responses", "analysis/13_virtual_perturbation/closed_loop/epas1_endothelial_program_response.png", 0)],
        cols=1,
    )
    compose(
        "Figure_5_drug_structure_closure",
        "Figure 5. Drug reversal and structure-guided prioritization",
        [
            ("A  DrugReflector reverse-signature ranking", "analysis/15_priority_figures/Figure_S12_DrugReflector_prioritization.png", 120),
            ("B  Five-conformer Vina affinities and co-crystal RMSD QC", "analysis/15_priority_figures/Figure_S16_EPAS1_main_docking_QC.png", 90),
            ("C  Representative ChimeraX EPAS1–ligand complexes", "analysis/15_priority_figures/Figure_S23_EPAS1_complex_binding_modes.png", 100),
        ],
        cols=1,
    )


if __name__ == "__main__":
    main()
