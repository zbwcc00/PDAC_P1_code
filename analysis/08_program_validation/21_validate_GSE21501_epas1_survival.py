from __future__ import annotations

import csv
import gzip
import math
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm

BASE = Path(r"D:/PDAC_P1")
RAW = BASE / "data/01_bulk/GSE21501/raw"
OUT = BASE / "analysis/08_program_validation/GSE21501"
OUT.mkdir(parents=True, exist_ok=True)

MATRIX = RAW / "GSE21501_series_matrix_full.txt.gz"
ANNOT = RAW / "GPL4133.annot.gz"
ENDOTHELIAL = [
    "EPAS1", "PECAM1", "VWF", "EMCN", "ESAM", "KDR", "ENG", "CDH5",
    "RAMP2", "CA4", "PLVAP", "RGCC", "SPARCL1", "MMRN2", "ANGPT2", "IL33",
]


def read_header_metadata(path: Path):
    samples = None
    rows = []
    with gzip.open(path, "rt", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                samples = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_characteristics_ch2"):
                rows.append(next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:])
            elif line.startswith("!series_matrix_table_begin"):
                break
    if samples is None:
        raise RuntimeError("GSE21501 sample identifiers were not found")
    metadata = pd.DataFrame({"sample": samples})
    # GEO occasionally concatenates multiple characteristics into one row
    # (e.g., tumor type for early samples and OS time for later samples).
    # Parse each cell independently so keys remain aligned to their sample.
    for values in rows:
        if len(values) != len(samples):
            continue
        for index, value in enumerate(values):
            if ":" not in value:
                continue
            key, parsed_value = value.split(":", 1)
            key = key.strip().lower()
            if not key:
                continue
            if key not in metadata:
                metadata[key] = ""
            metadata.loc[index, key] = parsed_value.strip()
    return metadata


def read_matrix(path: Path):
    lines = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            lines.append(line)
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, line in enumerate(lines[start:], start=start) if line.startswith("!series_matrix_table_end"))
    return pd.read_csv(StringIO("".join(lines[start:end])), sep="\t")


def read_annotation(path: Path):
    lines = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                break
        header = next(csv.reader([next(handle).rstrip("\n")], delimiter="\t"))
        for line in handle:
            if line.startswith("!platform_table_end"):
                break
            lines.append(line)
    ann = pd.read_csv(StringIO("".join(lines)), sep="\t", names=header,
                      header=None, dtype=str, low_memory=False)
    ann = ann[["ID", "Gene symbol"]].rename(columns={"ID": "ID_REF", "Gene symbol": "symbol"})
    ann["symbol"] = ann["symbol"].fillna("").str.split("///").str[0].str.strip()
    return ann


def cox_univariate(time, event, covariate):
    keep = np.isfinite(time) & np.isfinite(event) & np.isfinite(covariate) & (time > 0)
    time, event, covariate = time[keep], event[keep], covariate[keep]
    order = np.argsort(-time)
    time, event, covariate = time[order], event[order], covariate[order]

    def objective(beta):
        b = float(beta[0])
        eta = np.clip(b * covariate, -50, 50)
        risk = np.exp(eta)
        cumulative = np.cumsum(risk)
        loglik = np.sum(event * (eta - np.log(cumulative)))
        return -loglik

    result = minimize(objective, x0=np.array([0.0]), method="BFGS")
    beta = float(result.x[0])
    eta = np.clip(beta * covariate, -50, 50)
    risk = np.exp(eta)
    cumulative = np.cumsum(risk)
    cumulative_x = np.cumsum(risk * covariate)
    cumulative_x2 = np.cumsum(risk * covariate * covariate)
    variance = np.sum(event * (cumulative_x2 / cumulative - (cumulative_x / cumulative) ** 2))
    se = math.sqrt(1.0 / variance) if variance > 0 else float("nan")
    z = beta / se if np.isfinite(se) and se > 0 else float("nan")
    p = 2 * norm.sf(abs(z)) if np.isfinite(z) else float("nan")
    return {"n": int(keep.sum()), "events": int(event.sum()), "beta": beta,
            "HR": math.exp(beta), "SE": se, "P": p, "converged": bool(result.success)}


def km_curve(time, event, group):
    out = []
    at_risk = len(time)
    survival = 1.0
    for t in np.sort(np.unique(time[event == 1])):
        deaths = int(np.sum((time == t) & (event == 1)))
        censored = int(np.sum((time == t) & (event == 0)))
        out.append((t, survival))
        survival *= 1 - deaths / at_risk
        at_risk -= deaths + censored
    if not out:
        return np.array([]), np.array([])
    t, s = zip(*out)
    return np.asarray(t), np.asarray(s)


def main():
    metadata = read_header_metadata(MATRIX)
    expr = read_matrix(MATRIX)
    ann = read_annotation(ANNOT)
    expr["ID_REF"] = expr["ID_REF"].astype(str)
    ann["ID_REF"] = ann["ID_REF"].astype(str)
    expr = expr.merge(ann, on="ID_REF", how="left")
    expr["symbol"] = expr["symbol"].fillna("")
    expr = expr[expr["symbol"].isin(ENDOTHELIAL)].copy()
    sample_cols = [c for c in expr.columns if c in metadata["sample"].tolist()]
    gene = expr.groupby("symbol", as_index=True)[sample_cols].median()

    # Select samples with explicit overall-survival time and event metadata.
    clinical = metadata.copy()
    clinical["os_time"] = pd.to_numeric(clinical.get("os time", "").astype(str).str.extract(r"([0-9.]+)")[0], errors="coerce")
    clinical["os_event"] = pd.to_numeric(clinical.get("os event", "").astype(str).str.extract(r"([01])")[0], errors="coerce")
    clinical = clinical.dropna(subset=["os_time", "os_event"])
    clinical = clinical[clinical["os_time"] > 0].copy()
    clinical = clinical[clinical["sample"].isin(gene.columns)]
    clinical["os_event"] = clinical["os_event"].astype(int)
    clinical = clinical.drop_duplicates("sample")
    clinical = clinical.set_index("sample").loc[:, ["os_time", "os_event"]]

    available = [g for g in ENDOTHELIAL if g in gene.index]
    z = gene.loc[available, clinical.index].T
    z = (z - z.mean(axis=0)) / z.std(axis=0, ddof=1)
    clinical["EPAS1_expression"] = z["EPAS1"] if "EPAS1" in z else np.nan
    program_genes = [g for g in ENDOTHELIAL if g != "EPAS1" and g in z.columns]
    clinical["endothelial_program"] = z[program_genes].mean(axis=1) if program_genes else np.nan

    rows = []
    for feature in ["EPAS1_expression", "endothelial_program"]:
        result = cox_univariate(clinical["os_time"].to_numpy(float), clinical["os_event"].to_numpy(float), clinical[feature].to_numpy(float))
        result["feature"] = feature
        rows.append(result)
    results = pd.DataFrame(rows)[["feature", "n", "events", "HR", "SE", "P", "converged"]]
    results.to_csv(OUT / "GSE21501_EPAS1_external_survival.tsv", sep="\t", index=False)
    clinical.to_csv(OUT / "GSE21501_EPAS1_external_survival_input.tsv", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6), constrained_layout=True)
    for ax, feature, label in zip(axes, ["EPAS1_expression", "endothelial_program"], ["EPAS1", "Endothelial program"]):
        values = clinical[feature].to_numpy(float)
        cutoff = np.nanmedian(values)
        for group, color, name in [(values <= cutoff, "#0072B2", "Low"), (values > cutoff, "#D55E00", "High")]:
            t, s = km_curve(clinical.loc[group, "os_time"].to_numpy(float), clinical.loc[group, "os_event"].to_numpy(int), group)
            if len(t):
                ax.step(np.r_[0, t], np.r_[1, s], where="post", color=color, label=name)
        ax.set_title(label)
        ax.set_xlabel("Overall survival (months)")
        ax.set_ylabel("Kaplan–Meier survival")
        ax.set_ylim(0, 1.05)
        ax.legend(frameon=False, fontsize=8)
    fig.savefig(OUT / "GSE21501_EPAS1_external_survival_KM.png", dpi=300)
    fig.savefig(OUT / "GSE21501_EPAS1_external_survival_KM.pdf")

    report = [
        "# GSE21501 targeted external validation",
        "",
        f"- Expression platform: GPL4133 Agilent whole-human-genome array; {len(sample_cols)} matrix samples.",
        f"- Survival-eligible samples after explicit OS time/event parsing: n={len(clinical)}; events={int(clinical['os_event'].sum())}.",
        f"- EPAS1 probe/gene coverage: {'available' if 'EPAS1' in z.columns else 'not available'}; endothelial program genes available: {len(program_genes)}/{len(ENDOTHELIAL)-1}.",
        "- Analysis: prespecified EPAS1 expression and endothelial-program scores, z-standardized within this cohort; no model retraining or feature selection.",
        "- Cox estimates are univariate external associations and are not used to claim causal prognosis.",
        "",
        results.to_string(index=False),
        "",
        "Interpretation: this cohort is an independent clinical-stage/survival sensitivity validation. It does not validate spatial localization, drug efficacy, or EPAS1 causality.",
    ]
    (OUT / "GSE21501_EPAS1_external_validation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(OUT / "GSE21501_EPAS1_external_survival.tsv")


if __name__ == "__main__":
    main()
