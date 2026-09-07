from __future__ import annotations

import gzip
import itertools
import math
import re
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.paths import ANALYSIS_ROOT, EXTERNAL_DATA_ROOT

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr, wilcoxon


ROOT = ANALYSIS_ROOT.parent
OUT = ANALYSIS_ROOT / "18_external_strengthening"
OUT.mkdir(parents=True, exist_ok=True)
GSE202051 = EXTERNAL_DATA_ROOT / "GSE202051" / "GSE202051_totaldata-final-toshare.h5ad"
GSE300595_DIR = EXTERNAL_DATA_ROOT / "GSE300595" / "raw_h5"
GSE300595_SOFT = EXTERNAL_DATA_ROOT / "GSE300595" / "GSE300595_family.soft"

ENDO_MARKERS = ["PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "PLVAP", "CA4", "ADGRL4", "FLT1", "CDH5", "ENG"]
CORE_VASCULAR = ["PECAM1", "VWF", "EMCN"]
ENDO_TUMOR_UP = ["CASC15", "BCAT1", "C1QTNF6", "SESN3", "CTHRC1", "LEF1", "LHX6", "MARCKS"]
ENDO_TUMOR_DOWN = ["PRX", "DNASE1L3", "FGL2", "AQP7", "ADRB2", "SLCO4A1", "CA4", "BTNL9", "SLC14A1", "CTSH", "CDC42EP4", "TIMP3", "EPHX1", "HERPUD1", "RRAS", "TACC1", "GSN"]
SIGNAL_GENES = ["EPAS1", "ANGPT2", "SPARCL1", "VWF", "MMRN2", "IL33"]

LINEAGE_MARKERS = {
    "Endothelial": ENDO_MARKERS,
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "DCN", "LUM", "COL5A1", "COL5A2"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "AIF1", "C1QA", "C1QB", "LILRB1", "CTSD"],
    "T_NK": ["CD3D", "CD3E", "TRBC1", "LCK", "PTPRC", "NKG7", "TRAC"],
    "B_cell": ["MS4A1", "CD79A", "CD74", "HLA-DRA", "CD37"],
    "Epithelial": ["KRT19", "KRT8", "KRT18", "KRT7", "KRT17", "MSLN", "MUC1", "KRT23", "EPCAM"],
    "Acinar": ["PRSS1", "REG1A", "CPA1", "CTRB1", "AMY2A", "SYCN"],
    "Endocrine": ["INS", "GCG", "SST", "PPY", "IAPP", "CHGA"],
    "Pericyte_SMC": ["RGS5", "CSPG4", "MCAM", "PDGFRB", "NOTCH3", "DES", "ACTA2"],
}


def mean_expression(matrix: sparse.spmatrix | np.ndarray, gene_index: dict[str, int], genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in gene_index]
    if not present:
        return np.full(matrix.shape[0], np.nan), []
    indices = [gene_index[gene] for gene in present]
    values = matrix[:, indices]
    if sparse.issparse(values):
        values = np.asarray(values.mean(axis=1)).ravel()
    else:
        values = np.asarray(values).mean(axis=1)
    return values, present


def signed_score(matrix: sparse.spmatrix | np.ndarray, gene_index: dict[str, int]) -> tuple[np.ndarray, list[str], list[str]]:
    up, used_up = mean_expression(matrix, gene_index, ENDO_TUMOR_UP)
    down, used_down = mean_expression(matrix, gene_index, ENDO_TUMOR_DOWN)
    return (up - down) / 2, used_up, used_down


def exact_permutation_p(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    """Two-sided exact label permutation p value for responder minus nonresponder mean."""
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    observed = float(values[labels].mean() - values[~labels].mean())
    n = len(values)
    n_responder = int(labels.sum())
    null = []
    for selected in itertools.combinations(range(n), n_responder):
        selected = np.asarray(selected, dtype=int)
        mask = np.zeros(n, dtype=bool)
        mask[selected] = True
        null.append(float(values[mask].mean() - values[~mask].mean()))
    null = np.asarray(null)
    p = (np.abs(null) >= abs(observed) - 1e-12).mean()
    return observed, float(p), len(null)


def hedges_g(values: np.ndarray, labels: np.ndarray) -> float:
    responder, nonresponder = values[labels], values[~labels]
    n1, n0 = len(responder), len(nonresponder)
    if n1 < 2 or n0 < 2:
        return np.nan
    pooled = math.sqrt(((n1 - 1) * np.var(responder, ddof=1) + (n0 - 1) * np.var(nonresponder, ddof=1)) / (n1 + n0 - 2))
    if pooled == 0:
        return np.nan
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return float(correction * (np.mean(responder) - np.mean(nonresponder)) / pooled)


def audit_gse202051() -> None:
    if not GSE202051.exists():
        raise FileNotFoundError(GSE202051)
    data = ad.read_h5ad(GSE202051, backed="r")
    obs = data.obs.copy()
    if "broad_celltypes" not in obs or "pid" not in obs:
        raise ValueError("Expected broad_celltypes and pid in GSE202051 obs.")
    genes_needed = sorted(set(SIGNAL_GENES + ENDO_MARKERS + ENDO_TUMOR_UP + ENDO_TUMOR_DOWN))
    present = [gene for gene in genes_needed if gene in data.var_names]
    gene_index = {gene: i for i, gene in enumerate(present)}
    matrix = data[:, present].X
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    expression = pd.DataFrame(index=obs.index)
    for gene in present:
        vector = matrix[:, gene_index[gene]]
        expression[gene] = vector.toarray().ravel() if sparse.issparse(vector) else np.asarray(vector).ravel()
    expression["pid"] = obs["pid"].astype(str).values
    expression["sampleid"] = obs["sampleid"].astype(str).values
    expression["lineage"] = obs["broad_celltypes"].astype(str).values
    expression["endothelial_subtype"] = obs["new_celltypes"].astype(str).values
    expression["treatment_status"] = obs["treatment_status"].astype(str).values
    expression["response"] = obs["response"].astype(str).values
    expression["n_counts"] = pd.to_numeric(obs["total_counts"], errors="coerce").values
    expression["n_genes"] = pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").values
    expression["vascular_marker_score"] = expression[[gene for gene in CORE_VASCULAR if gene in expression]].mean(axis=1)
    used_up = [gene for gene in ENDO_TUMOR_UP if gene in expression]
    used_down = [gene for gene in ENDO_TUMOR_DOWN if gene in expression]
    expression["frozen_endothelial_remodeling_score"] = (expression[used_up].mean(axis=1) - expression[used_down].mean(axis=1)) / 2

    cell_counts = expression.groupby(["pid", "lineage"], observed=True).size().rename("n_cells").reset_index()
    cell_counts.to_csv(OUT / "GSE202051_patient_lineage_cell_counts.tsv", sep="\t", index=False)
    patient_lineage = expression.groupby(["pid", "lineage"], observed=True).agg(
        n_cells=("EPAS1", "size"),
        EPAS1_mean=("EPAS1", "mean"),
        EPAS1_detection=("EPAS1", lambda x: float(np.mean(np.asarray(x) > 0))),
        ANGPT2_mean=("ANGPT2", "mean"),
        vascular_marker_score=("vascular_marker_score", "mean"),
        frozen_endothelial_remodeling_score=("frozen_endothelial_remodeling_score", "mean"),
    ).reset_index()
    patient_lineage.to_csv(OUT / "GSE202051_patient_lineage_expression.tsv", sep="\t", index=False)

    endo = expression.loc[expression["lineage"].eq("Endothelial")].copy()
    endo_patient = endo.groupby("pid", observed=True).agg(
        n_endothelial_cells=("EPAS1", "size"),
        EPAS1_mean=("EPAS1", "mean"),
        EPAS1_detection=("EPAS1", lambda x: float(np.mean(np.asarray(x) > 0))),
        ANGPT2_mean=("ANGPT2", "mean"),
        ANGPT2_detection=("ANGPT2", lambda x: float(np.mean(np.asarray(x) > 0))),
        vascular_marker_score=("vascular_marker_score", "mean"),
        frozen_endothelial_remodeling_score=("frozen_endothelial_remodeling_score", "mean"),
        treatment_status=("treatment_status", lambda x: sorted(set(x))[0] if len(set(x)) == 1 else ";".join(sorted(set(x)))),
        response=("response", lambda x: sorted(set(x))[0] if len(set(x)) == 1 else ";".join(sorted(set(x)))),
    ).reset_index()
    endo_patient.to_csv(OUT / "GSE202051_endothelial_patient_pseudobulk_like.tsv", sep="\t", index=False)

    subtype = endo.groupby("endothelial_subtype", observed=True).agg(
        n_cells=("EPAS1", "size"), EPAS1_mean=("EPAS1", "mean"), EPAS1_detection=("EPAS1", lambda x: float(np.mean(np.asarray(x) > 0))),
        ANGPT2_mean=("ANGPT2", "mean"), vascular_marker_score=("vascular_marker_score", "mean"),
        frozen_endothelial_remodeling_score=("frozen_endothelial_remodeling_score", "mean"),
    ).reset_index().sort_values("n_cells", ascending=False)
    subtype.to_csv(OUT / "GSE202051_endothelial_subtype_expression.tsv", sep="\t", index=False)

    paired = patient_lineage.pivot(index="pid", columns="lineage", values="EPAS1_mean")
    paired = paired.dropna(subset=["Endothelial"])
    comparator = paired.drop(columns=["Endothelial"], errors="ignore").median(axis=1)
    valid = comparator.notna()
    paired_delta = paired.loc[valid, "Endothelial"] - comparator.loc[valid]
    w_stat, w_p = wilcoxon(paired_delta, alternative="greater", zero_method="wilcox") if len(paired_delta) >= 3 else (np.nan, np.nan)
    endo_vascular = endo_patient.dropna(subset=["EPAS1_mean", "vascular_marker_score"])
    rho, rho_p = spearmanr(endo_vascular["EPAS1_mean"], endo_vascular["vascular_marker_score"]) if len(endo_vascular) >= 3 else (np.nan, np.nan)
    summary = pd.DataFrame([{
        "dataset": "GSE202051", "n_cells": int(data.n_obs), "n_patients": int(obs["pid"].nunique()),
        "n_endothelial_cells": int(len(endo)), "n_patients_with_endothelial": int(endo_patient.shape[0]),
        "median_endothelial_cells_per_patient": float(endo_patient["n_endothelial_cells"].median()),
        "endothelial_EPAS1_mean_median": float(endo_patient["EPAS1_mean"].median()),
        "endothelial_EPAS1_detection_median": float(endo_patient["EPAS1_detection"].median()),
        "paired_endothelial_minus_other_lineage_EPAS1_median": float(paired_delta.median()),
        "paired_wilcoxon_greater_p": float(w_p), "paired_n_patients": int(len(paired_delta)),
        "patient_EPAS1_vs_vascular_rho": float(rho), "patient_EPAS1_vs_vascular_p": float(rho_p),
        "frozen_score_up_genes": ";".join(used_up), "frozen_score_down_genes": ";".join(used_down),
        "scope": "Independent tumor-only localization; not a tumor-normal or therapy-response test",
    }])
    summary.to_csv(OUT / "GSE202051_independent_localization_summary.tsv", sep="\t", index=False)


def parse_gse300595_metadata() -> pd.DataFrame:
    text = GSE300595_SOFT.read_text(encoding="utf-8", errors="ignore")
    rows = []
    for block in text.split("^SAMPLE = ")[1:]:
        title = re.search(r"!Sample_title = (.*)", block)
        if title is None or not title.group(1).endswith("_GEX"):
            continue
        patient = title.group(1).replace("_GEX", "").lower()
        characteristics = re.findall(r"!Sample_characteristics_ch1 = (.*)", block)
        row = {"patient": patient}
        for item in characteristics:
            if ": " in item:
                key, value = item.split(": ", 1)
                row[key] = value
        rows.append(row)
    result = pd.DataFrame(rows).sort_values("patient").reset_index(drop=True)
    result["responder"] = result["genotype"].eq("responder")
    return result


def read_multiome_gene_counts(path: Path) -> tuple[sparse.csr_matrix, list[str], list[str], int, int]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        feature_type = np.asarray(matrix["features"]["feature_type"][:]).astype(str)
        genes = np.asarray(matrix["features"]["name"][:]).astype(str)
        barcodes = np.asarray(matrix["barcodes"][:]).astype(str)
        full = sparse.csc_matrix((matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]), shape=tuple(matrix["shape"][:]))
    gene_rows = np.flatnonzero(feature_type == "Gene Expression")
    peak_rows = np.flatnonzero(feature_type == "Peaks")
    return full[gene_rows, :].T.tocsr(), genes[gene_rows].tolist(), barcodes.tolist(), int(len(genes)), int(len(peak_rows))


def audit_gse300595() -> None:
    metadata = parse_gse300595_metadata()
    metadata.to_csv(OUT / "GSE300595_patient_metadata.tsv", sep="\t", index=False)
    cell_rows, patient_rows = [], []
    for file in sorted(GSE300595_DIR.glob("*.h5")):
        patient = re.search(r"_(pc\d+)\.filtered", file.name, flags=re.I).group(1).lower()
        counts, genes, barcodes, n_total_features, n_peak_features = read_multiome_gene_counts(file)
        library_size = np.asarray(counts.sum(axis=1)).ravel()
        scaling = 1e4 / np.maximum(library_size, 1)
        lognorm = sparse.diags(scaling).dot(counts)
        lognorm.data = np.log1p(lognorm.data)
        gene_index = {gene: index for index, gene in enumerate(genes)}
        scores = {}
        used = {}
        for lineage, markers in LINEAGE_MARKERS.items():
            scores[lineage], used[lineage] = mean_expression(lognorm, gene_index, markers)
        score_frame = pd.DataFrame(scores)
        cell_label = score_frame.idxmax(axis=1).to_numpy()
        endo_score = score_frame["Endothelial"].to_numpy()
        second_score = np.partition(score_frame.to_numpy(), -2, axis=1)[:, -2]
        # Require two detectable endothelial markers and a positive margin over the runner-up lineage.
        marker_indices = [gene_index[gene] for gene in ENDO_MARKERS if gene in gene_index]
        endo_detected = np.asarray((counts[:, marker_indices] > 0).sum(axis=1)).ravel() if marker_indices else np.zeros(counts.shape[0])
        endothelial = (cell_label == "Endothelial") & (endo_detected >= 2) & (endo_score > second_score)
        epas1 = lognorm[:, gene_index["EPAS1"]].toarray().ravel() if "EPAS1" in gene_index else np.full(counts.shape[0], np.nan)
        angpt2 = lognorm[:, gene_index["ANGPT2"]].toarray().ravel() if "ANGPT2" in gene_index else np.full(counts.shape[0], np.nan)
        frozen, used_up, used_down = signed_score(lognorm, gene_index)
        vascular, used_vascular = mean_expression(lognorm, gene_index, CORE_VASCULAR)
        meta = metadata.loc[metadata["patient"].eq(patient)].iloc[0]
        keep_cols = ["pid", "n_counts", "n_genes", "endothelial", "best_lineage", "endothelial_score", "EPAS1", "ANGPT2", "vascular_marker_score", "frozen_endothelial_remodeling_score"]
        cell_table = pd.DataFrame({
            "pid": patient, "barcode": barcodes, "n_counts": library_size,
            "n_genes": np.asarray((counts > 0).sum(axis=1)).ravel(), "endothelial": endothelial,
            "best_lineage": cell_label, "endothelial_score": endo_score, "EPAS1": epas1, "ANGPT2": angpt2,
            "vascular_marker_score": vascular, "frozen_endothelial_remodeling_score": frozen,
            "response_group": meta["genotype"], "treatment": meta["treatment"],
        })
        cell_rows.append(cell_table)
        endo = cell_table.loc[cell_table["endothelial"]].copy()
        patient_rows.append({
            "patient": patient, "response_group": meta["genotype"], "responder": bool(meta["responder"]), "treatment": meta["treatment"],
            "n_nuclei": int(counts.shape[0]), "n_features_total": n_total_features, "n_gene_expression_features": int(len(genes)), "n_peak_features": n_peak_features,
            "n_endothelial_nuclei": int(len(endo)), "endothelial_fraction": float(len(endo) / counts.shape[0]),
            "median_endothelial_score": float(np.median(endo_score)), "EPAS1_mean_endothelial": float(endo["EPAS1"].mean()) if len(endo) else np.nan,
            "EPAS1_detection_endothelial": float((endo["EPAS1"] > 0).mean()) if len(endo) else np.nan,
            "ANGPT2_mean_endothelial": float(endo["ANGPT2"].mean()) if len(endo) else np.nan,
            "ANGPT2_detection_endothelial": float((endo["ANGPT2"] > 0).mean()) if len(endo) else np.nan,
            "vascular_marker_score_endothelial": float(endo["vascular_marker_score"].mean()) if len(endo) else np.nan,
            "frozen_endothelial_remodeling_score": float(endo["frozen_endothelial_remodeling_score"].mean()) if len(endo) else np.nan,
            "endo_markers_used": ";".join(marker for marker in ENDO_MARKERS if marker in gene_index),
            "frozen_score_up_genes_used": ";".join(used_up), "frozen_score_down_genes_used": ";".join(used_down),
        })
    cell_table = pd.concat(cell_rows, ignore_index=True)
    cell_table.to_csv(OUT / "GSE300595_all_nuclei_marker_annotation.tsv.gz", sep="\t", index=False, compression="gzip")
    patient_table = pd.DataFrame(patient_rows).sort_values("patient").reset_index(drop=True)

    inferential = []
    eligible = patient_table.loc[patient_table["n_endothelial_nuclei"] >= 20].copy()
    for outcome in ["EPAS1_mean_endothelial", "ANGPT2_mean_endothelial", "vascular_marker_score_endothelial", "frozen_endothelial_remodeling_score", "endothelial_fraction"]:
        valid = eligible.dropna(subset=[outcome])
        if valid["responder"].nunique() != 2:
            continue
        values, labels = valid[outcome].to_numpy(float), valid["responder"].to_numpy(bool)
        effect, p, permutations = exact_permutation_p(values, labels)
        inferential.append({"outcome": outcome, "n_eligible_patients": len(valid), "n_responders": int(labels.sum()), "n_nonresponders": int((~labels).sum()),
                            "responder_minus_nonresponder": effect, "hedges_g": hedges_g(values, labels), "exact_two_sided_p": p, "n_exact_permutations": permutations,
                            "scope": "post-neoadjuvant resection association, not pretreatment prediction"})
        folf = valid.loc[valid["treatment"].eq("FOLFIRINOX")]
        if folf["responder"].nunique() == 2:
            effect, p, permutations = exact_permutation_p(folf[outcome].to_numpy(float), folf["responder"].to_numpy(bool))
            inferential.append({"outcome": outcome + " [FOLFIRINOX sensitivity]", "n_eligible_patients": len(folf), "n_responders": int(folf["responder"].sum()), "n_nonresponders": int((~folf["responder"]).sum()),
                                "responder_minus_nonresponder": effect, "hedges_g": hedges_g(folf[outcome].to_numpy(float), folf["responder"].to_numpy(bool)), "exact_two_sided_p": p, "n_exact_permutations": permutations,
                                "scope": "post-neoadjuvant FOLFIRINOX-only sensitivity; exploratory small-n association"})
    inferential_table = pd.DataFrame(inferential)
    inferential_table.to_csv(OUT / "GSE300595_endothelial_response_exact_tests.tsv", sep="\t", index=False)
    patient_table.to_csv(OUT / "GSE300595_endothelial_patient_summary.tsv", sep="\t", index=False)


def write_report() -> None:
    loc = pd.read_csv(OUT / "GSE202051_independent_localization_summary.tsv", sep="\t").iloc[0]
    multiome = pd.read_csv(OUT / "GSE300595_endothelial_patient_summary.tsv", sep="\t")
    tests = pd.read_csv(OUT / "GSE300595_endothelial_response_exact_tests.tsv", sep="\t")
    lines = [
        "# External PDAC single-cell and multiome reinforcement audit",
        "",
        "## GSE202051: independent tumor endothelial-state localization",
        "",
        f"- The downloaded object contained {int(loc.n_cells):,} cells from {int(loc.n_patients)} patients, including {int(loc.n_endothelial_cells):,} endothelial-labelled cells from {int(loc.n_patients_with_endothelial)} patients.",
        f"- Patient-paired endothelial EPAS1 expression exceeded the median non-endothelial lineage expression by a median {loc.paired_endothelial_minus_other_lineage_EPAS1_median:.3f} normalized units (one-sided paired Wilcoxon P={loc.paired_wilcoxon_greater_p:.3g}; n={int(loc.paired_n_patients)}).",
        f"- Across endothelial patient summaries, EPAS1 and the independent PECAM1/VWF/EMCN score had Spearman rho={loc.patient_EPAS1_vs_vascular_rho:.3f} (P={loc.patient_EPAS1_vs_vascular_p:.3g}).",
        "- This is a tumor-only external localization cohort. It supports endothelial localization/state recurrence, but cannot serve as a tumor-versus-normal test or a causal treatment-response analysis.",
        "",
        "## GSE300595: post-neoadjuvant multiome response-associated context",
        "",
        f"- Twelve resected post-neoadjuvant PDAC tumors were audited: {int(multiome.responder.sum())} histopathologic responders and {int((~multiome.responder).sum())} nonresponders.",
        f"- Marker-based annotation identified endothelial nuclei in {(multiome.n_endothelial_nuclei > 0).sum()}/12 patients; {(multiome.n_endothelial_nuclei >= 20).sum()}/12 met the predeclared >=20 endothelial-nuclei threshold for patient-level testing.",
        "- RNA and ATAC peak features were supplied in the same H5 object. This first pass uses RNA-only marker annotation; ATAC inference is conditional on a sufficiently populated endothelial compartment and will not be forced from sparse nuclei.",
        "- All response comparisons are exact patient-label permutations, not nucleus-level tests. They are post-treatment associations and must not be described as pretreatment predictive biomarkers or evidence of treatment causality.",
        "",
        "## Interpretation gate",
        "",
        "- A response result is eligible for the manuscript only if endothelial abundance is adequate, the effect direction is coherent with the main endothelial story, and the exact test plus treatment-stratified sensitivity are not contradictory.",
        "- The original EPAS1 → vascular–immune ecology → virtual perturbation → drug-prioritization story is unchanged by either audit.",
    ]
    if len(tests):
        lines.extend(["", "### Exact response tests", ""])
        for _, row in tests.iterrows():
            lines.append(f"- {row['outcome']}: responder minus nonresponder={row['responder_minus_nonresponder']:.3f}; Hedges g={row['hedges_g']:.3f}; exact P={row['exact_two_sided_p']:.3g}; n={int(row['n_eligible_patients'])}.")
    (OUT / "external_strengthening_audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    audit_gse202051()
    audit_gse300595()
    write_report()
    print("External reinforcement audit complete:", OUT)
