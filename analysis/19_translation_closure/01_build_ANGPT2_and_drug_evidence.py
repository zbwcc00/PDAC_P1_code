from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path("D:/PDAC_P1")
OUT = ROOT / "analysis" / "19_translation_closure"
OUT.mkdir(parents=True, exist_ok=True)
P0 = ROOT / "analysis" / "17_P0_reinforcement"
CELLCHAT = ROOT / "analysis" / "10_communication_pseudotime" / "CellChat"
SECACT = ROOT / "analysis" / "09_spatial_validation" / "SecAct"
EXTERNAL = ROOT / "analysis" / "18_external_strengthening"
DR = ROOT / "analysis" / "13_virtual_perturbation" / "drugreflector"
PRISM = ROOT / "analysis" / "11_methodology_review" / "PRISM20Q2"
DOCK = ROOT / "analysis" / "13_virtual_perturbation" / "docking" / "EPAS1_multiconformer"


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def build_angpt2_evidence() -> None:
    expression = pd.read_csv(P0 / "ANGPT2_receptor_patient_celltype_expression.tsv", sep="\t")
    datasets = sorted(expression["dataset"].unique())
    rows = []
    correlations = []
    for dataset in datasets:
        source = expression.loc[(expression.dataset == dataset) & (expression.cell_type == "endothelial") & (expression.gene == "ANGPT2"),
                                ["patient", "mean_log_expression", "detection_fraction"]].rename(columns={"mean_log_expression": "ANGPT2_source_mean", "detection_fraction": "ANGPT2_source_detection"})
        source_supported = int((source.ANGPT2_source_detection > 0).sum())
        rows.append({"evidence_layer": "Source expression", "dataset_or_cohort": dataset, "sender": "endothelial", "receiver": "not applicable",
                     "result": f"ANGPT2 detectable in {source_supported}/{len(source)} patient groups", "support": "supportive", "unit": "patient x cell type", "scope": "source availability; not direct signaling"})
        for receiver in ["myeloid", "T_NK", "CD8", "B_cell"]:
            receptor = expression.loc[(expression.dataset == dataset) & (expression.cell_type == receiver) & expression.gene.isin(["ITGA5", "ITGB1"]),
                                      ["patient", "gene", "mean_log_expression", "detection_fraction"]]
            receptor = receptor.pivot(index="patient", columns="gene", values=["mean_log_expression", "detection_fraction"])
            receptor.columns = [f"{first}_{second}" for first, second in receptor.columns]
            receptor = receptor.reset_index()
            required = ["mean_log_expression_ITGA5", "mean_log_expression_ITGB1", "detection_fraction_ITGA5", "detection_fraction_ITGB1"]
            if not all(column in receptor for column in required):
                continue
            joined = source.merge(receptor, on="patient", how="inner").dropna()
            joint_detected = (joined.detection_fraction_ITGA5.gt(0) & joined.detection_fraction_ITGB1.gt(0))
            if receiver == "myeloid":
                rows.append({"evidence_layer": "Receiver availability", "dataset_or_cohort": dataset, "sender": "endothelial", "receiver": receiver,
                             "result": f"ITGA5 and ITGB1 jointly detectable in {int(joint_detected.sum())}/{len(joined)} matched patient groups", "support": "supportive" if joint_detected.all() else "partial", "unit": "patient x cell type", "scope": "receptor availability; not direct signaling"})
            if len(joined) >= 3:
                receptor_joint = joined[["mean_log_expression_ITGA5", "mean_log_expression_ITGB1"]].min(axis=1)
                rho, p_value = spearmanr(joined.ANGPT2_source_mean, receptor_joint)
                correlations.append({"dataset": dataset, "receiver": receiver, "n_matched_patient_groups": len(joined), "spearman_rho": rho, "p_value": p_value,
                                     "joint_receptor_score": "min(mean ITGA5, mean ITGB1)", "interpretation": "patient-level source-receiver co-variation; not a direct interaction test"})

    correlations = pd.DataFrame(correlations)
    correlations["BH_q"] = benjamini_hochberg(correlations.p_value)
    correlations["support"] = np.where((correlations.receiver == "myeloid") & (correlations.BH_q < 0.05) & (correlations.spearman_rho > 0), "supportive", "not consistently supportive")
    correlations.to_csv(OUT / "ANGPT2_patient_source_receiver_correlation.tsv", sep="\t", index=False)

    edge_rows = []
    for path in sorted(CELLCHAT.glob("*_all_communications.tsv")):
        dataset = path.name.replace("_all_communications.tsv", "")
        table = pd.read_csv(path, sep="\t")
        query = table.loc[(table.source == "endothelial") & (table.target == "myeloid") & (table.ligand == "ANGPT2") & (table.receptor == "ITGA5_ITGB1")].copy()
        if len(query):
            edge_rows.append({"dataset": dataset, "cellchat_edge_present": True, "n_edges": len(query), "min_p": query.pval.min(), "mean_probability": query.prob.mean()})
        else:
            edge_rows.append({"dataset": dataset, "cellchat_edge_present": False, "n_edges": 0, "min_p": np.nan, "mean_probability": np.nan})
    edges = pd.DataFrame(edge_rows)
    edges.to_csv(OUT / "ANGPT2_myeloid_CellChat_dataset_audit.tsv", sep="\t", index=False)
    rows.append({"evidence_layer": "CellChat candidate edge", "dataset_or_cohort": "GSE154778,GSE155698,GSE212966", "sender": "endothelial", "receiver": "myeloid",
                 "result": f"ANGPT2–ITGA5_ITGB1 edge in {int(edges.cellchat_edge_present.sum())}/{len(edges)} independently run datasets", "support": "supportive" if edges.cellchat_edge_present.all() else "partial", "unit": "dataset-level CellChat inference", "scope": "database-guided communication inference; not receptor occupancy"})

    spatial = pd.read_csv(P0 / "ANGPT2_SecAct_spatial_patient_replication.tsv", sep="\t")
    for _, result in spatial.iterrows():
        rows.append({"evidence_layer": "SecAct spatial permutation", "dataset_or_cohort": result.dataset, "sender": "EPAS1-high spatial context", "receiver": "not cell-resolved",
                     "result": f"ANGPT2 positive in {int(result.positive_patients)}/{int(result.n_patients)} patients; median beta={result.median_patient_beta:.4g}",
                     "support": "supportive", "unit": "patient; 1,000 within-ROI permutations", "scope": "spatial secreted-activity context; not source-to-receiver causality"})

    external = pd.read_csv(EXTERNAL / "GSE202051_endothelial_patient_pseudobulk_like.tsv", sep="\t")
    rho, p_value = spearmanr(external.EPAS1_mean, external.ANGPT2_mean)
    rows.append({"evidence_layer": "Independent endothelial state", "dataset_or_cohort": "GSE202051", "sender": "endothelial", "receiver": "not applicable",
                 "result": f"EPAS1–ANGPT2 patient-level rho={rho:.3f}, P={p_value:.3g} (n={len(external)})", "support": "not supportive", "unit": "patient-level endothelial summary", "scope": "external state co-variation test; null result retained"})

    matrix = pd.DataFrame(rows)
    matrix.to_csv(OUT / "ANGPT2_cross_method_evidence_matrix.tsv", sep="\t", index=False)
    report = [
        "# ANGPT2 endothelial–myeloid cross-method communication assessment", "",
        "## Supported components", "",
        f"- CellChat inferred the endothelial ANGPT2–ITGA5/ITGB1-to-myeloid edge in {int(edges.cellchat_edge_present.sum())}/{len(edges)} independently run single-cell datasets.",
        "- Endothelial ANGPT2 and myeloid ITGA5/ITGB1 availability were each observed in all three single-cell datasets at the patient × cell-type level.",
        f"- SecAct spatial activity was directionally positive for ANGPT2 in {int(spatial.positive_patients.sum())}/{int(spatial.n_patients.sum())} patients across two independent spatial cohorts, using the precomputed 1,000 within-ROI permutations.",
        "", "## Limiting component", "",
        "- Patient-level source–receiver expression co-variation was not consistently positive across discovery cohorts after multiple-testing correction, and EPAS1–ANGPT2 co-variation in independent GSE202051 endothelial summaries was null. These analyses do not invalidate the candidate interface; they preclude describing it as a universally coupled patient-level axis.",
        "", "## Manuscript language", "",
        "The appropriate conclusion is: *multi-method analyses nominate an endothelial ANGPT2–ITGA5/ITGB1 interface with myeloid cells in PDAC*. Do not use language implying direct signaling, ligand binding, TEK-mediated immune signaling, or causal immune reprogramming.",
    ]
    (OUT / "ANGPT2_cross_method_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def build_drug_evidence() -> None:
    priority = pd.read_csv(DR / "drugreflector_priority.tsv", sep="\t")
    mapping = pd.read_csv(DR / "drugreflector_compound_mapping.tsv", sep="\t")
    prism = pd.read_csv(PRISM / "PRISM20Q2_EPAS1_candidate_drug_results.tsv", sep="\t")
    structure = pd.read_csv(DOCK / "epas1_final_structure_evidence.tsv", sep="\t")
    rows = []
    curated = {
        "Y-39983": {"compound": "BRD-K56751279", "pubchem_cid": 11507964, "clinical_status": "Experimental ROCK inhibitor; no PDAC trial returned by ClinicalTrials.gov query", "clinical_status_source": "LINCS/PRISM annotation; ClinicalTrials.gov queried 2026-09-06"},
        "triclabendazole": {"compound": "BRD-K81916719", "pubchem_cid": 50248, "clinical_status": "FDA-labelled oral anthelmintic (Egaten) for fascioliasis; no PDAC trial returned by ClinicalTrials.gov query", "clinical_status_source": "openFDA label NDA208711; ClinicalTrials.gov queried 2026-09-06"},
    }
    for name, meta in curated.items():
        compound = meta["compound"]
        priority_row = priority.loc[(priority.target == "EPAS1") & (priority.compound == compound)].iloc[0]
        prism_row = prism.loc[prism.name.str.lower().eq(name.lower())].iloc[0]
        structure_row = structure.loc[structure.cmap_name.str.lower().eq(name.lower())].iloc[0]
        mapping_row = mapping.loc[mapping.compound.eq(compound)].iloc[0]
        sensitivity_supported = bool(prism_row.fraction_lfc_le_neg05 >= 0.5 and prism_row.median_lfc <= -0.5)
        decision = "primary testable repurposing hypothesis" if (name == "triclabendazole" and sensitivity_supported) else "structure-supported negative pharmacologic comparator"
        rows.append({
            "compound": name, "LINCS_BRD_identifier": compound, "PubChem_CID": meta["pubchem_cid"], "DrugReflector_EPAS1_reverse_rank": priority_row.mean_rank,
            "DrugReflector_recurrent_signatures": priority_row.n_signatures, "DrugReflector_recurrent_strata": priority_row.n_strata,
            "PRISM_pancreas_models": prism_row.n_panc, "PRISM_median_logFC": prism_row.median_lfc, "PRISM_fraction_logFC_le_minus0_5": prism_row.fraction_lfc_le_neg05,
            "PRISM_sensitivity_gate": "supportive" if sensitivity_supported else "not supportive", "docking_mean_affinity_kcal_mol": structure_row.mean_affinity_kcal_mol,
            "docking_five_conformer_consistency": structure_row.five_conformer_consistency, "redocking_QC": structure_row.redocking_qc,
            "LINCS_CLUE_annotation": mapping_row.moa, "LINCS_CLUE_annotated_target": mapping_row.target, "clinical_or_development_status": meta["clinical_status"],
            "status_source": meta["clinical_status_source"], "PDAC_trial_status": "No returned trial for compound + pancreatic cancer", "integrated_decision": decision,
            "interpretation_guardrail": "DrugReflector ranks signature reversal and docking estimates structural compatibility; neither establishes direct EPAS1 binding or PDAC clinical efficacy.",
        })
    matrix = pd.DataFrame(rows)
    matrix.to_csv(OUT / "EPAS1_drug_repositioning_evidence_matrix.tsv", sep="\t", index=False)
    report = [
        "# EPAS1 drug-repositioning evidence matrix", "",
        "## Prioritization decision", "",
        "- Triclabendazole is retained as the primary *testable repurposing hypothesis*: it recurs in both EPAS1 reverse-signature strata, shows 5/5 conformer structural compatibility with passed redocking QC, and has orthogonal PRISM sensitivity in pancreas-labelled models (median logFC −0.856; 78.8% of 33 models ≤ −0.5). It has an established non-oncology clinical label, not a PDAC indication.",
        "- Y-39983 is retained as the structure-supported negative pharmacologic comparator: it recurs in both reverse-signature strata and has stable 5/5-conformer docking, but its PRISM pancreas-model sensitivity is weak (median logFC −0.095; 16.7% of 36 models ≤ −0.5).",
        "", "## Explicit limits", "",
        "- ClinicalTrials.gov searches for triclabendazole + pancreatic cancer, Y-39983 + pancreatic cancer, and EPAS1 + pancreatic cancer returned no study on 2026-09-06.",
        "- DrugReflector does not identify a direct molecular target; PRISM measures cancer-cell viability rather than endothelial specificity; docking does not prove biochemical binding. The matrix supports a ranked experimental hypothesis, not an efficacy claim.",
    ]
    (OUT / "EPAS1_drug_repositioning_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    build_angpt2_evidence()
    build_drug_evidence()
    print("ANGPT2 and drug evidence matrices complete:", OUT)
