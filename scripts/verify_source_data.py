#!/usr/bin/env python3
"""Check that deposited source-data tables contain the manuscript's key values."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data"


def load(name: str) -> pd.DataFrame:
    path = SOURCE / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def close(actual: float, expected: float, tolerance: float = 1e-8) -> bool:
    return abs(float(actual) - expected) <= tolerance * max(1.0, abs(expected))


def main() -> int:
    checks: list[tuple[str, bool]] = []

    scale = load("Figure1_dataset_scale.tsv")
    checks.append(("discovery cohort sample count", int(scale["samples"].sum()) == 66))
    checks.append(("discovery cohort cell count", int(scale["total_cells"].sum()) == 224546))

    independent = load("Figure2_GSE202051_localization_summary.tsv").iloc[0]
    checks.append(("GSE202051 patient count", int(independent["n_patients"]) == 43))
    checks.append(("GSE202051 endothelial cell count", int(independent["n_endothelial_cells"]) == 19258))
    checks.append(("GSE202051 paired Wilcoxon P", close(independent["paired_wilcoxon_greater_p"], 1.136868e-13, 1e-6)))
    checks.append(("GSE202051 EPAS1-vascular rho", close(independent["patient_EPAS1_vs_vascular_rho"], 0.7783147085472667)))

    conditional = load("Table2_conditional_model_results.tsv")
    endothelial = conditional.loc[conditional["outcome"] == "endothelial_neighbor_score"].iloc[0]
    tnk = conditional.loc[conditional["outcome"] == "T_NK_neighbor_score"].iloc[0]
    myeloid = conditional.loc[conditional["outcome"] == "myeloid_neighbor_score"].iloc[0]
    checks.append(("GSE282302 endothelial conditional beta", close(endothelial["median_conditional_beta"], 0.09109849, 1e-5)))
    checks.append(("GSE282302 T/NK conditional q", close(tnk["BH_q_across_spatial_endpoints"], 0.0258789, 1e-5)))
    checks.append(("GSE282302 myeloid conditional q", close(myeloid["BH_q_across_spatial_endpoints"], 0.5652669, 1e-5)))

    pharmacology = load("Figure5_pharmacology_source_data.tsv")
    tric = pharmacology.loc[pharmacology["name"] == "triclabendazole"].iloc[0]
    checks.append(("triclabendazole PRISM model count", int(tric["n_panc"]) == 33))
    checks.append(("triclabendazole PRISM median logFC", close(tric["median_lfc"], -0.856462467267, 1e-6)))

    docking = load("Figure5_docking_source_data.tsv")
    checks.append(("docking local paths removed", not {"pose_file", "log_file"}.intersection(docking.columns)))

    failed = [label for label, passed in checks if not passed]
    for label, passed in checks:
        print(f"{'PASS' if passed else 'FAIL'}\t{label}")
    if failed:
        print(f"{len(failed)} source-data checks failed", file=sys.stderr)
        return 1
    print(f"All {len(checks)} source-data checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
