from pathlib import Path
from PIL import Image, ImageOps, ImageDraw

ROOT = Path(r"D:/第二篇大论文/analysis/15_priority_figures")
OUT = ROOT / "audit_2026-09-05"

groups = {
    "supp_S6_S14": sorted(ROOT.glob("Figure_S[6-9]_*.png")) + sorted(ROOT.glob("Figure_S1[0-4]_*.png")),
    "supp_S15_S23": sorted(ROOT.glob("Figure_S15_*.png")) + sorted(ROOT.glob("Figure_S16_*.png")) + sorted(ROOT.glob("Figure_S17_*.png")) + sorted(ROOT.glob("Figure_S18_*.png")) + sorted(ROOT.glob("Figure_S19_*.png")) + sorted(ROOT.glob("Figure_S20_*.png")) + sorted(ROOT.glob("Figure_S21_*.png")) + sorted(ROOT.glob("Figure_S22_*.png")) + sorted(ROOT.glob("Figure_S23_*.png")),
    "supp_S24_S25": sorted((ROOT / "figure2_supplement").glob("Figure_S24*.png")) + sorted((ROOT / "figure2_supplement").glob("Figure_S25*.png")),
    "main": sorted((ROOT / "main_figures_unified").glob("Figure_*.png")),
}

for name, files in groups.items():
    files = [f for f in files if f.exists()]
    if not files:
        continue
    thumbs = []
    for f in files:
        im = Image.open(f).convert("RGB")
        im.thumbnail((700, 500))
        canvas = Image.new("RGB", (720, 550), "white")
        x = (720 - im.width) // 2
        canvas.paste(im, (x, 30))
        d = ImageDraw.Draw(canvas)
        d.text((10, 8), f.name, fill="black")
        thumbs.append(canvas)
    cols = 2
    rows = (len(thumbs) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * 720, rows * 550), "#dddddd")
    for i, im in enumerate(thumbs):
        sheet.paste(im, ((i % cols) * 720, (i // cols) * 550))
    sheet.save(OUT / f"{name}_contact.png", dpi=(150, 150))
