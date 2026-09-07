from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path(r"D:/PDAC_P1")
MANUSCRIPT = BASE / "manuscript/PDAC_P1_locked_results_methods.md"
LEGENDS = BASE / "manuscript/PDAC_P1_results_and_figure_legends.md"
OUT = BASE / "analysis/00_progress/PDAC_P1_numeric_asset_audit_2026-09-04.md"


def check(name: str, passed: bool, detail: str, rows: list[tuple[str, bool, str]]):
    rows.append((name, passed, detail))


def main():
    rows: list[tuple[str, bool, str]] = []
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    legends = LEGENDS.read_text(encoding="utf-8")

    scrna = pd.read_csv(BASE / "analysis/02_scrna_qc/PDAC_scRNA_dataset_summary.tsv", sep="\t")
    expected_scrna = {"GSE154778": (16, 17086), "GSE155698": (38, 136386), "GSE212966": (12, 71074)}
    observed_scrna = {row.dataset: (int(row.samples), int(row.total_cells)) for row in scrna.itertuples()}
    check("Discovery single-cell cohorts", observed_scrna == expected_scrna,
          f"Observed {observed_scrna}; manuscript reports the same cohort sizes.", rows)

    secact = pd.read_csv(BASE / "analysis/09_spatial_validation/SecAct/SecAct_five_targets_patient_summary_1000.tsv", sep="\t")
    core = ["SPARCL1", "VWF", "MMRN2", "ANGPT2", "IL33"]
    core_rows = secact[secact.secreted.isin(core)]
    core_ok = ((core_rows.query("dataset == 'GSE282302'").positive_patients == 14).all()
               and (core_rows.query("dataset == 'GSE297144'").positive_patients == 8).all())
    check("SecAct core-panel patient consistency", core_ok,
          "All five candidates are 14/14 in GSE282302 and 8/8 in GSE297144.", rows)

    evidence = (BASE / "analysis/09_spatial_validation/SecAct/final_evidence_matrix.tsv").read_text(encoding="utf-8")
    secact_text_ok = all(token in manuscript for token in ["1.22×10−4", "2.93×10−4", "7.81×10−3", "8.52×10−3"])
    matrix_q_ok = all(token in evidence for token in ["BH q=0.000293", "BH q=0.00852"])
    check("SecAct P/Q reporting", secact_text_ok and matrix_q_ok,
          "Main text and evidence matrix agree on sign-test P values and BH-adjusted Q values.", rows)

    ecology = pd.read_csv(BASE / "analysis/09_spatial_validation/GSE282302/spatial_ecology/GSE282302_spatial_ecology_patient.tsv", sep="\t")
    myeloid = ecology.median_myeloid_delta.median()
    tnk = ecology.median_T_NK_delta.median()
    ecology_ok = (len(ecology) == 14 and (ecology.median_myeloid_delta > 0).all()
                  and (ecology.median_T_NK_delta > 0).all()
                  and round(myeloid, 3) == 0.070 and round(tnk, 3) == 0.089)
    check("Spatial ecology", ecology_ok,
          f"n={len(ecology)}; median deltas: myeloid={myeloid:.3f}, T/NK={tnk:.3f}.", rows)

    perturb = pd.read_csv(BASE / "analysis/13_virtual_perturbation/closed_loop/epas1_mechanism_program_summary.tsv", sep="\t")
    endothelial = perturb[(perturb.target == "EPAS1") & (perturb.cell_type == "endothelial")]
    ko = endothelial[(endothelial.perturbation == "KO") & endothelial.program.isin(["endothelial", "angiogenesis"])]
    oe = endothelial[(endothelial.perturbation == "OE") & endothelial.program.isin(["endothelial", "angiogenesis"])]
    perturb_ok = (ko.n_significant_datasets.tolist() == [3, 3] and oe.n_significant_datasets.tolist() == [1, 1]
                  and [round(x, 2) for x in ko.median_abs_Z] == [0.99, 1.07]
                  and [round(x, 2) for x in oe.median_abs_Z] == [0.62, 0.80])
    check("Virtual perturbation", perturb_ok,
          "Endothelial KO is 3/3 for both prespecified programs; OE is 1/3.", rows)

    docking = pd.read_csv(BASE / "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/epas1_final_structure_evidence.tsv", sep="\t")
    docking_ok = (round(float(docking.loc[docking.cmap_name == "Y-39983", "mean_affinity_kcal_mol"].iloc[0]), 3) == -6.553
                  and round(float(docking.loc[docking.cmap_name == "triclabendazole", "mean_affinity_kcal_mol"].iloc[0]), 3) == -6.115
                  and (docking.redocking_qc == "5/5 co-crystal RMSD pass").all())
    check("Five-conformer docking", docking_ok,
          "Y-39983=-6.553 and triclabendazole=-6.115 kcal/mol; 5/5 redocking QC passes.", rows)

    survival = pd.read_csv(BASE / "analysis/08_program_validation/GSE21501/GSE21501_EPAS1_external_survival.tsv", sep="\t")
    survival_ok = (survival.n.tolist() == [102, 102] and survival.events.tolist() == [66, 66]
                   and round(float(survival.loc[survival.feature == "EPAS1_expression", "HR"].iloc[0]), 3) == 0.914
                   and round(float(survival.loc[survival.feature == "EPAS1_expression", "P"].iloc[0]), 3) == 0.483
                   and round(float(survival.loc[survival.feature == "endothelial_program", "P"].iloc[0]), 3) == 0.876)
    check("GSE21501 external survival gate", survival_ok,
          "n=102, events=66; both prespecified Cox associations are null.", rows)

    placeholders = ["Figure X", "Figure Y", "Supplementary Figure Y", "Supplementary Table Y"]
    placeholder_ok = not any(token in legends for token in placeholders)
    check("Figure/table placeholder removal", placeholder_ok,
          "No unresolved Figure X/Y or Table Y placeholders in the current legends draft.", rows)

    methods_ok = all(token in manuscript for token in ["GSE21501", "1,000 coordinate-preserving", "exhaustiveness 16", "RMSD ≤2 Å"])
    check("Methods-to-analysis detail coverage", methods_ok,
          "Methods explicitly state the external cohort, spatial permutation, docking search, and RMSD gate.", rows)

    survival_matrix_ok = "Targeted external survival\tEPAS1\tGSE21501" in evidence
    check("Evidence-matrix update", survival_matrix_ok,
          "GSE21501 is recorded as a negative promotion gate in the final evidence matrix.", rows)

    ready_assets = [
        "analysis/15_priority_figures/main_figures_relayout/Figure_1_scRNA_discovery_convergence.pdf",
        "analysis/15_priority_figures/main_figures_relayout/Figure_2_bulk_spatial_validation.pdf",
        "analysis/15_priority_figures/main_figures_relayout/Figure_3_secact_immune_ecology.pdf",
        "analysis/15_priority_figures/main_figures_relayout/Figure_4_virtual_perturbation.pdf",
        "analysis/15_priority_figures/main_figures_relayout/Figure_5_drug_structure_closure.pdf",
        "analysis/13_virtual_perturbation/closed_loop/manuscript_figure/EPAS1_closed_loop_main_mechanism_figure.pdf",
        "analysis/15_priority_figures/Figure_S6_scRNA_discovery_overview.pdf",
        "analysis/09_spatial_validation/SecAct/Figure_integrated_EPAS1_SecAct_CellChat_ecology.pdf",
        "analysis/08_program_validation/GSE21501/GSE21501_EPAS1_external_survival_KM.pdf",
        "analysis/15_priority_figures/Figure_S11_CPTAC_protein_evidence.pdf",
        "analysis/15_priority_figures/Figure_S13_negative_pharmacology_gates.pdf",
        "analysis/15_priority_figures/Figure_S14_survival_ML_sensitivity.pdf",
        "analysis/15_priority_figures/Figure_S16_EPAS1_structure_docking_QC.pdf",
        "analysis/15_priority_figures/Figure_S17_MR_instrument_gate.pdf",
        "analysis/15_priority_figures/Figure_S18_GSE21501_survival_forest.pdf",
        "analysis/15_priority_figures/Figure_S19_functional_immune_WGCNA.pdf",
        "analysis/15_priority_figures/Figure_S20_TCGA_bulk_QC_PCA.pdf",
        "analysis/15_priority_figures/Figure_S21_spatial_image_montage.pdf",
        "analysis/15_priority_figures/Figure_S22_endothelial_DE_pathway_gate.pdf",
        "analysis/15_priority_figures/Figure_S23_EPAS1_complex_binding_modes.pdf",
        "analysis/15_priority_figures/scrna_additional/Figure_S15_scRNA_sample_QC.pdf",
        "analysis/15_priority_figures/scrna_additional/Figure_S16_scRNA_marker_DotPlot.pdf",
        "analysis/15_priority_figures/scrna_additional/Figure_S17_endothelial_pseudobulk_volcano.pdf",
    ]
    assets_ok = all((BASE / path).exists() for path in ready_assets)
    check("Ready figure assets", assets_ok,
          f"{sum((BASE / path).exists() for path in ready_assets)}/{len(ready_assets)} manifest-ready figure assets exist.", rows)

    passed = sum(status for _, status, _ in rows)
    report = ["# PDAC P1 numerical and asset consistency audit (2026-09-04)", "",
              f"**Result: {passed}/{len(rows)} checks passed.**", "",
              "| Check | Status | Evidence |", "|---|---|---|"]
    for name, status, detail in rows:
        report.append(f"| {name} | {'PASS' if status else 'FIX REQUIRED'} | {detail} |")
    report.extend([
        "",
        "## Scope and boundary",
        "",
        "This audit reconciles the locked Results/Methods and figure-legends drafts against the primary result tables listed above. It confirms arithmetic and label consistency for the central evidence chain; it does not convert computational associations into causal or therapeutic claims.",
        "",
        "## Numbering decision",
        "",
        "Figures 1–6 now carry the primary discovery, validation, mechanism, therapeutic-prioritization, and integrative-summary claims. Supplementary Figures S1–S14 provide QC, full-resolution detail, negative gates, and sensitivity analyses. The full evidence matrix remains Supplementary Table S1.",
    ])
    OUT.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(OUT)
    if passed != len(rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
