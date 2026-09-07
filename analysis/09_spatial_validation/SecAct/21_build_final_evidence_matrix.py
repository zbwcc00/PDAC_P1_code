from pathlib import Path

import pandas as pd

BASE = Path(r"D:/PDAC_P1")
SEC = BASE / "analysis/09_spatial_validation/SecAct"
OUT = SEC / "final_evidence_matrix.tsv"
REPORT = SEC / "final_evidence_matrix.md"

rows = []

def add(layer, target, cohort, direction, consistency, perm_fdr, cellchat,
        ecology, grade, interpretation, source):
    rows.append({
        "evidence_layer": layer,
        "target": target,
        "cohort_or_dataset": cohort,
        "direction": direction,
        "patient_level_consistency": consistency,
        "permutation_or_FDR": perm_fdr,
        "CellChat_support": cellchat,
        "spatial_ecology_support": ecology,
        "evidence_grade": grade,
        "interpretation": interpretation,
        "source_file": source,
    })

summary = pd.read_csv(SEC / "SecAct_five_targets_patient_summary_1000.tsv", sep="\t")

def bh_adjust(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    adjusted = [1.0] * len(values)
    running = 1.0
    n = len(values)
    for rank, idx in reversed(list(enumerate(order, start=1))):
        running = min(running, values[idx] * n / rank)
        adjusted[idx] = running
    return adjusted

summary["exact_sign_q"] = bh_adjust(summary["exact_sign_p"].astype(float).tolist())
for _, r in summary.iterrows():
    target = r["secreted"]
    cohort = r["dataset"]
    positive = f"{int(r['positive_patients'])}/{int(r['n_patients'])}"
    p = float(r["exact_sign_p"])
    q = float(r["exact_sign_q"])
    if target == "IGFBP7":
        grade = "C"
        interp = "directional sensitivity signal; not retained in the core spatial panel"
    else:
        grade = "A"
        interp = "cross-cohort patient-level spatial signal"
    add(
        "SecAct spatial permutation", target, cohort, "positive",
        positive, f"exact sign p={p:.3g}; BH q={q:.3g}; 1000 coordinate-preserving permutations/ROI",
        "ANGPT2 only" if target == "ANGPT2" else "none detected",
        "EPAS1-linked context (not target-specific)", grade, interp,
        "SecAct_five_targets_patient_summary_1000.tsv",
    )

add("CellChat", "ANGPT2", "GSE282302", "endothelial source -> immune receivers",
    "3 receiver classes; 5 total edges", "permutation p not applicable",
    "ITGA5_ITGB1: myeloid, CD8, B cell", "supports myeloid/T/NK neighborhood context",
    "B", "candidate-specific communication support; association, not causal signaling",
    "SecAct_candidate_CellChat_edges.tsv")

add("Spatial ecology", "EPAS1", "GSE282302", "EPAS1-high > EPAS1-low",
    "14/14 patients for myeloid and T/NK programs", "patient sign test reported in spatial ecology output",
    "indirect", "median delta myeloid=0.070; T/NK=0.089",
    "B", "immune neighborhood context for the EPAS1 axis; CAF association unstable",
    "SecAct_spatial_ecology_context_summary.tsv")

add("MR primary", "EPAS1", "GTEx v8 -> PDAC GWAS", "not estimable",
    "0 genome-wide GTEx instruments", "failed instrument/outcome-overlap gate",
    "not applicable", "not applicable", "D",
    "no causal MR claim; report as an explicit negative gate", "MR_instrument_gate_summary.tsv")
add("MR orthogonal", "TACC1; HERPUD1", "BLUEPRINT/Lepik -> PDAC GWAS", "exploratory only",
    "cis signals present but tissue/outcome mismatch", "no primary valid instrument set",
    "not applicable", "not applicable", "C",
    "sensitivity evidence only; not used to establish causality", "MR_instrument_gate_summary.tsv")

add("Protein abundance", "TACC1; MARCKS; HERPUD1", "CPTAC-PAAD", "detected",
    "140, 140, and 93 samples non-missing", "not applicable",
    "not applicable", "orthogonal to spatial context", "B",
    "protein-layer detectability supports biological assayability", "CPTAC_protein_detection_summary.tsv")
add("Protein abundance", "EPAS1", "HPA", "protein evidence; tissue/cell specificity not detected",
    "evidence at protein level", "not applicable", "not applicable", "not applicable", "C",
    "protein existence is supported, but localization evidence is weak", "HPA_protein_layer_summary.tsv")
add("Clinical protein survival", "TACC1", "CPTAC-PAAD", "higher protein associated with outcome",
    "n=134; 76 events", "continuous Cox p=0.018; median split p=0.334",
    "not applicable", "not applicable", "B",
    "continuous association is supportive but requires external validation", "CPTAC_protein_survival_cox.tsv")

add("DepMap CRISPR", "EPAS1; MARCKS; HERPUD1; TACC1", "47 PDAC cell models",
    "weak dependency", "strong-dependency fraction=0 for all four", "not applicable",
    "not applicable", "not applicable", "D", "does not support a robust tumor-cell autonomous dependency",
    "DepMap_PDAC_candidate_dependency_summary.tsv")
add("PRISM pharmacology", "Y-39983; triclabendazole", "PDAC cell-line subset",
    "association not significant", "n=33-36; all FDR > 0.70", "not applicable",
    "not applicable", "not applicable", "D", "in-vitro drug-response association is not reproducible enough for promotion",
    "DepMap_PRISM_association_results.tsv")
add("Drug-response ML", "Y-39983; triclabendazole", "PRISM + molecular features",
    "no predictive gain", "100x repeated 5-fold CV; median R2 < 0; retention gate failed",
    "outer-CV gate failed", "not applicable", "not applicable", "D", "SHAP/ML retained as transparent negative sensitivity analysis",
    "drug_ml_repeated_CV_summary.tsv")

add("Virtual perturbation", "EPAS1", "GSE154778/GSE155698/GSE212966",
    "KO/OE alters endothelial/angiogenesis programs", "endothelial program significant in 3/3 KO datasets",
    "multi-dataset program FDR gate", "not applicable", "not applicable", "B",
    "multi-dataset computational perturbation supports mechanistic prioritization",
    "epas1_mechanism_program_summary.tsv")
add("DrugReflector", "EPAS1", "LINCS signatures; composite/endothelial strata",
    "reverse signature recurrence", "Y-39983 and triclabendazole recurrent in 2/2 strata",
    "recurrent in 2/2 strata", "not applicable", "not applicable", "B", "prioritization evidence for hypothesis generation",
    "drugreflector_priority.tsv")
add("Molecular docking", "EPAS1", "5 PAS-B conformers",
    "favorable predicted binding", "Y-39983 and triclabendazole stable across 5/5 conformers",
    "5/5 co-crystal RMSD QC pass", "not applicable", "not applicable", "B", "structure-supported computational evidence; not direct binding validation",
    "epas1_final_structure_evidence.tsv")

df = pd.DataFrame(rows)
df.to_csv(OUT, sep="\t", index=False)

core = ["ANGPT2", "SPARCL1", "VWF", "MMRN2", "IL33"]
lines = [
    "# Final evidence matrix",
    "",
    "Core spatial candidates: **ANGPT2, SPARCL1, VWF, MMRN2, IL33**.",
    "IGFBP7 is retained as a sensitivity candidate because replication is incomplete in GSE297144 (5/8 positive patients).",
    "",
    "Evidence grades: A = replicated patient-level spatial evidence; B = orthogonal/supportive evidence; C = limited or sensitivity evidence; D = negative/failed gate.",
    "",
    f"Rows: {len(df)}; core candidates: {', '.join(core)}.",
    "",
    "The matrix is descriptive and does not convert computational associations into experimental causality.",
]
REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"wrote {OUT}")
print(f"wrote {REPORT}")
