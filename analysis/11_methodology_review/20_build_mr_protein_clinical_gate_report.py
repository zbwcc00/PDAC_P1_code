from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd

BASE = Path(r"D:\PDAC_P1")
MR = BASE / "analysis" / "11_methodology_review"
OUT = MR / "supplemental_gate"
OUT.mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def main():
    gtex = read_tsv(MR / "GTEx_v8_targeted" / "GTEx_PDAC_GWAS_harmonization_audit.tsv")
    blueprint = read_tsv(MR / "eqtl_catalogue_targeted" / "BLUEPRINT_targeted_cis_eqtl_summary.tsv")
    lepik = read_tsv(MR / "eqtl_catalogue_targeted" / "LEPIK_targeted_cis_eqtl_summary.tsv")
    rows = []
    for gene in ["EPAS1", "TACC1", "HERPUD1"]:
        g = gtex[gtex["gene"] == gene]
        b = blueprint[blueprint["target_gene"] == gene].iloc[0]
        l = lepik[lepik["target_gene"] == gene].iloc[0]
        rows.append({
            "gene": gene,
            "GTEx_overlap_n": len(g),
            "GTEx_overlap_min_eqtl_p": float(g["eqtl_p"].min()) if len(g) else None,
            "GTEx_genomewide_instrument": bool(len(g) and (g["eqtl_p"] < 5e-8).any()),
            "BLUEPRINT_cis_n": int(b["n"]),
            "BLUEPRINT_n_p5e8": int(b["n_p5e8"]),
            "BLUEPRINT_n_p1e5": int(b["n_p1e5"]),
            "LEPIK_cis_n": int(l["n"]),
            "LEPIK_n_p5e8": int(l["n_p5e8"]),
            "LEPIK_n_p1e5": int(l["n_p1e5"]),
            "LEPIK_max_F": float(l["max_F"]),
            "primary_MR_status": "fail" if len(g) == 0 or not (g["eqtl_p"] < 5e-8).any() else "conditional",
            "orthogonal_MR_status": "exploratory_only: tissue_mismatch_or_no_outcome_overlap",
        })
    gate = pd.DataFrame(rows)
    gate.to_csv(OUT / "MR_instrument_gate_summary.tsv", sep="\t", index=False)

    hpa = read_tsv(MR / "HPA" / "HPA_target_summary.tsv")
    hpa.to_csv(OUT / "HPA_protein_layer_summary.tsv", sep="\t", index=False)

    report = [
        "# Supplemental evidence gate: MR, protein, clinical and DepMap/PRISM",
        "",
        "## MR",
        "- The prespecified primary instrument gate requires cis-eQTL P<5×10⁻⁸, allele harmonization, F>10 and LD clumping.",
        "- GTEx–PDAC coordinate overlaps for EPAS1/TACC1/HERPUD1 do not provide a genome-wide-significant cis-eQTL instrument; no primary MR estimate is therefore reported.",
        "- Lepik provides strong whole-blood cis-eQTLs for all three genes, but tissue mismatch and missing/limited outcome overlap mean these remain orthogonal exploratory sensitivity data, not endothelial causal evidence.",
        "- BLUEPRINT has no genome-wide-significant cis-eQTL for these genes and no usable PDAC outcome overlap.",
        "",
        "## Protein layer",
        "- HPA gene-level records were downloaded for EPAS1, TACC1, MARCKS and HERPUD1.",
        "- EPAS1 is annotated as evidence at protein level but is not detected in the HPA summary's protein cell-type/tissue distribution; this does not support a strong pancreatic endothelial protein claim.",
        "- HERPUD1 shows cell-type enrichment/group-enriched tissue distribution; MARCKS is broadly detected with tissue/cell-type enhancement; TACC1 is broadly detected with low specificity.",
        "- No CPTAC pancreatic quantitative protein table is currently present; therefore HPA is supportive annotation, not PDAC-specific protein validation.",
        "",
        "## Clinical external validation",
        "- GSE62165 contains 131 HG-U219 samples with tissue and grouped-stage labels. EPAS1 expression and the endothelial score showed no grouped-stage association (Kruskal–Wallis P=0.424 and P=0.953, respectively).",
        "- The GSE62165 series matrix has no survival time/censoring endpoint. GSE21501 is a mixed tissue/reference-RNA expression series without a usable survival endpoint in the downloaded matrix.",
        "- These cohorts can support expression/stage sensitivity analyses, but not a prognostic model or survival claim.",
        "",
        "## DepMap/PRISM",
        "- DepMap 24Q4 CRISPRGeneEffect is available and includes 47 PDAC models after lineage filtering.",
        "- EPAS1, MARCKS and HERPUD1 have median gene effects −0.037, −0.071 and −0.073, with 0/47 models at the strong-dependency threshold ≤−1; this is weak descriptive dependency evidence.",
        "- PRISM compound-sensitivity data are not present in the local D-drive data directory; no PRISM efficacy claim is made. DrugReflector remains a computational ranking only.",
        "",
        "## Decision",
        "- Keep EPAS1 as the primary perturbational/structural target, but do not label it MR-proven, protein-validated or PDAC-selectively dependent.",
        "- Treat MR as a formally failed/insufficient evidence layer for the current datasets; report this transparently rather than adding weak instruments.",
        "- HPA, GSE62165 and DepMap are useful boundary-setting support; the main causal gap remains unresolved without stronger cis-pQTL/eQTL, PDAC protein data, or experimental validation.",
    ]
    (OUT / "supplemental_evidence_gate_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("supplemental evidence gate report complete")


if __name__ == "__main__":
    main()
