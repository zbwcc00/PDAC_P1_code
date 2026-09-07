from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import importlib.util
spec = importlib.util.spec_from_file_location("epas_validation", HERE.parent / "20_epas1_external_spatial_validation.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

def main():
    scores = pd.read_csv(mod.ROOT / "GSE327056_all_spot_program_scores_global.tsv", sep="\t")
    scaling = pd.read_csv(mod.ROOT / "GSE327056_program_gene_global_scaling.tsv", sep="\t").set_index("gene")
    program = mod.read_program_genes()
    down_n = max(1, len(program["down"]))
    sample_info = {"A1": ("GSM9647219", "adjacent_tumor"), "B1": ("GSM9647220", "tumor"), "C1": ("GSM9647221", "tumor_stroma"), "D1": ("GSM9647222", "normal_pancreas")}
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), constrained_layout=True)
    for j, (tag, (gsm, context)) in enumerate(sample_info.items()):
        pos, epas = mod.read_visium(tag, gsm)
        one = scores[scores["tag"] == tag].set_index("barcode").loc[pos["barcode"]].reset_index()
        epas_z = (epas - float(scaling.loc["EPAS1", "global_mean"])) / float(scaling.loc["EPAS1", "global_sd"])
        one["endothelial_score_exEPAS1"] = one["endothelial_score"] + epas_z / down_n
        x, y = pos["pxl_col_in_fullres"].to_numpy(), pos["pxl_row_in_fullres"].to_numpy()
        for ax, values, title, cmap in [(axes[0, j], epas, f"{tag} {context}\nEPAS1", "viridis"), (axes[1, j], one["endothelial_score_exEPAS1"].to_numpy(), "Endothelial program (EPAS1 excluded)", "magma")]:
            order = np.argsort(values)
            sc = ax.scatter(x[order], y[order], c=values[order], s=1.2, cmap=cmap, linewidths=0)
            ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_title(title, fontsize=9); ax.axis("off")
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.01)
    fig.suptitle("EPAS1 spatial localization across PDAC sections", fontsize=14)
    fig.savefig(HERE / "EPAS1_spatial_localization_maps.png", dpi=350)
    fig.savefig(HERE / "EPAS1_spatial_localization_maps.pdf")
    print("maps_written")

if __name__ == "__main__":
    main()
