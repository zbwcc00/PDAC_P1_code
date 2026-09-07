from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASE = Path(r"D:\PDAC_P1")
ROOT = BASE / "analysis" / "13_virtual_perturbation"
DRUG = ROOT / "drugreflector"
DOCK = ROOT / "docking" / "EPAS1_multiconformer"
OUT = ROOT / "closed_loop"
OUT.mkdir(parents=True, exist_ok=True)


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    vp = read(ROOT / "virtual_perturbation_network_program_summary.tsv")
    epas = [r for r in vp if r["target"] == "EPAS1"]
    focus_programs = ["endothelial", "angiogenesis", "TNF_NFkB", "myeloid_inflammation", "CAF_ECM", "cytotoxic_TNK", "T_exhaustion"]
    mechanism = []
    for cell_type in sorted({r["cell_type"] for r in epas}):
        for perturbation in ["KO", "OE"]:
            for program in focus_programs:
                rows = [r for r in epas if r["cell_type"] == cell_type and r["perturbation"] == perturbation and r["program"] == program]
                if not rows:
                    continue
                sig = [int(r["n_sig_fdr"]) > 0 for r in rows]
                mechanism.append({
                    "target": "EPAS1", "cell_type": cell_type, "perturbation": perturbation, "program": program,
                    "n_datasets": len(rows), "n_significant_datasets": sum(sig),
                    "median_abs_Z": f"{np.median([float(r['mean_abs_Z']) for r in rows]):.4f}",
                    "median_mean_FC": f"{np.median([float(r['mean_FC']) for r in rows]):.6g}",
                    "evidence": "multi-dataset program response",
                })
    write(OUT / "epas1_mechanism_program_summary.tsv", mechanism)

    final = read(DOCK / "final_exhaustiveness16" / "final_exhaustiveness16_summary.tsv")
    mapping = {r["compound"]: r for r in read(DRUG / "drugreflector_compound_mapping_pubchem.tsv")}
    candidates = []
    for compound in ["BRD-K56751279", "BRD-K81916719"]:
        rows = [r for r in final if r["compound"] == compound]
        affinities = np.array([float(r["affinity_best_kcal_mol"]) for r in rows])
        m = mapping[compound]
        candidates.append({
            "target": "EPAS1", "compound": compound, "cmap_name": m.get("cmap_name", ""), "pubchem_cid": m.get("pubchem_cid", ""),
            "n_conformers": len(rows), "mean_affinity_kcal_mol": f"{affinities.mean():.4f}", "sd_affinity_kcal_mol": f"{affinities.std(ddof=1):.4f}",
            "min_affinity_kcal_mol": f"{affinities.min():.4f}", "max_affinity_kcal_mol": f"{affinities.max():.4f}",
            "drugreflector_evidence": "EPAS1 reverse top50 recurrent across composite/endothelial signatures",
            "structure_evidence": "five-conformer EPAS1 PAS-B docking; positive co-crystal redocking QC passed",
            "decision": "primary" if compound == "BRD-K56751279" else "secondary",
        })
    write(OUT / "epas1_closed_loop_candidates.tsv", candidates)

    # Main mechanism heatmap: endothelial and angiogenesis response under EPAS1 KO/OE.
    plot_rows = [r for r in mechanism if r["cell_type"] == "endothelial" and r["program"] in ["endothelial", "angiogenesis", "TNF_NFkB", "cytotoxic_TNK"]]
    labels = [f"{r['perturbation']}:{r['program']}" for r in plot_rows]
    values = [float(r["median_abs_Z"]) for r in plot_rows]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    colors = ["#b2182b" if r["perturbation"] == "KO" else "#2166ac" for r in plot_rows]
    ax.barh(labels, values, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Median |Z| (across datasets)")
    ax.set_title("EPAS1 perturbation response in endothelial programs")
    fig.tight_layout(); fig.savefig(OUT / "epas1_endothelial_program_response.png", dpi=300); plt.close(fig)

    # Final docking stability plot.
    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    names = ["Y-39983", "Triclabendazole"]
    values = [[float(r["affinity_best_kcal_mol"]) for r in final if r["compound"] == c] for c in ["BRD-K56751279", "BRD-K81916719"]]
    ax.boxplot(values, labels=names, showmeans=True)
    ax.set_ylabel("Vina affinity (kcal/mol)")
    ax.set_title("EPAS1 PAS-B five-conformer final docking")
    fig.tight_layout(); fig.savefig(OUT / "epas1_final_docking_stability.png", dpi=300); plt.close(fig)
    print("closed-loop integration complete")


if __name__ == "__main__":
    main()
