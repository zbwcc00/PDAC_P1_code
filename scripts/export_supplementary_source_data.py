#!/usr/bin/env python3
"""Export non-identifying result tables supporting Supplementary Figures S1-S32."""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import sys
from pathlib import Path

import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from config.paths import ANALYSIS_ROOT


OUTPUT_ROOT = REPOSITORY_ROOT / "source_data" / "supplementary"
PATH_COLUMN = re.compile(r"(?:^|_)(?:path|file|directory|dir|log)(?:$|_)", re.IGNORECASE)
ABSOLUTE_PATH = re.compile(r"(?:(?:[A-Za-z]:\\)|(?:/Users/)|(?:/home/))")

# filename, upstream analysis-relative path, manuscript panels, biological row unit, description
EXPORTS = [
    ("S01_sample_qc.tsv", "02_scrna_qc/PDAC_scRNA_sample_QC.tsv", "S1,S16", "GEO sample", "Sample-level single-cell quality-control metrics."),
    ("S01_patient_lineage_composition.tsv", "03_scrna_annotation/PDAC_scRNA_patient_lineage_summary.tsv", "S1,S15,S17", "patient-lineage summary", "Annotated discovery-cohort lineage composition."),
    ("S02_secact_evidence.tsv", "09_spatial_validation/SecAct/final_evidence_matrix.tsv", "S2,S22,S29", "candidate evidence row", "Spatial SecAct inference summary."),
    ("S02_cellchat_endothelial_consensus.tsv", "10_communication_pseudotime/CellChat_endothelial_consensus.tsv", "S2,S22,S29", "inferred cell-cell edge", "CellChat endothelial consensus; computational inference."),
    ("S03_S09_GSE21501_survival.tsv", "08_program_validation/GSE21501/GSE21501_EPAS1_external_survival.tsv", "S3,S9", "survival model", "External survival sensitivity estimates."),
    ("S04_cptac_protein_detection.tsv", "11_methodology_review/CPTAC_PAAD/CPTAC_protein_detection_summary.tsv", "S4", "protein target", "CPTAC protein-detection summary."),
    ("S04_cptac_protein_survival.tsv", "11_methodology_review/CPTAC_PAAD/CPTAC_protein_survival_cox.tsv", "S4", "survival model", "CPTAC survival sensitivity estimates."),
    ("S05_depmap_prism_association.tsv", "11_methodology_review/DepMap_PRISM_pharmacology/DepMap_CRISPR_PRISM_association_results.tsv", "S5", "gene-drug association", "DepMap/PRISM prioritization boundary analysis."),
    ("S05_drug_ml_cv_summary.tsv", "14_prediction/DepMap_PRISM_drug_ML/drug_ml_repeated_CV_summary.tsv", "S5", "model summary", "Repeated cross-validation summary for drug ML."),
    ("S06_survival_ml_performance.tsv", "14_prediction/epas1_survival_ml/survival_model_performance.tsv", "S6", "model summary", "TCGA/CPTAC survival-model performance summary."),
    ("S07_S24_docking_summary.tsv", "13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_exhaustiveness16_summary.tsv", "S7,S24", "structure-compound summary", "Five-conformer docking and redocking quality-control summary."),
    ("S07_S24_docking_contacts.tsv", "13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_contact_by_conformer.tsv", "S7,S14,S24", "structure-compound-contact", "Geometry-based docking contact candidates."),
    ("S08_mr_instrument_gate.tsv", "11_methodology_review/supplemental_gate/MR_instrument_gate_summary.tsv", "S8", "instrument assessment", "MR instrument availability and promotion gate."),
    ("S10_scRNA_immune_correlations.tsv", "12_functional_immune/scRNA_endothelial_immune_correlations_by_dataset.tsv", "S10", "cohort association", "Discovery-cohort immune correlations."),
    ("S10_hallmark_gsea.tsv", "12_functional_immune/TCGA_program_correlated_Hallmark_GSEA.tsv", "S10", "gene set", "TCGA Hallmark GSEA results."),
    ("S10_wgcna_module_traits.tsv", "12_functional_immune/WGCNA/TCGA_WGCNA_module_trait_correlations.tsv", "S10", "module-trait association", "TCGA WGCNA module-trait correlations."),
    ("S11_tcga_qc.tsv", "01_tcga_qc/TCGA-PAAD_primary_tumor_qc.tsv", "S11", "TCGA sample", "TCGA primary-tumor quality-control values."),
    ("S12_S26_GSE282302_roi_anchor.tsv", "09_spatial_validation/GSE282302/GSE282302_ROI_spatial_anchor_results.tsv", "S12,S21,S26", "ROI", "GSE282302 image-resolved spatial anchor results."),
    ("S13_S18_de_concordance.tsv", "06_pseudobulk_DE/celltype_DE_concordance.tsv", "S13,S18,S19", "cell-type/cohort result", "Cross-method pseudobulk differential-expression concordance."),
    ("S13_S19_target_prioritization.tsv", "08_program_validation/endothelial_target_prioritization.tsv", "S13,S19", "candidate target", "Endothelial target prioritization summary."),
    ("S13_pathway_gate.tsv", "15_priority_figures/S22_endothelial_DE_GO_KEGG_recurrent.tsv", "S13", "pathway", "Recurrent endothelial DE pathway gate."),
    ("S15_copykat_sample_summary.tsv", "04_cnv_copykat/PDAC_copykat_all_sample_summary.tsv", "S15", "sample/cell-type summary", "Exploratory CopyKAT/CNV sample summary."),
    ("S15_S17_lineage_score_summary.tsv", "03_scrna_annotation/PDAC_scRNA_patient_lineage_restricted_ecology.tsv", "S15,S17", "patient-lineage score", "Patient-level annotated lineage-score summary supporting annotation displays."),
    ("S20_GSE62452_candidate_results.tsv", "08_program_validation/GSE62452/GSE62452_paired_candidate_gene_results.tsv", "S20,S25", "gene contrast", "Paired GSE62452 candidate-gene validation."),
    ("S20_GSE71729_candidate_results.tsv", "08_program_validation/GSE71729/GSE71729_primary_vs_normal_candidate_gene_results.tsv", "S20,S25", "gene contrast", "GSE71729 primary-versus-normal candidate-gene validation."),
    ("S21_spatial_ecology_patient.tsv", "09_spatial_validation/GSE282302/spatial_ecology/GSE282302_spatial_ecology_patient.tsv", "S21", "patient", "Patient-level unadjusted spatial ecology results."),
    ("S22_ANGPT2_cellchat_edges.tsv", "17_P0_reinforcement/ANGPT2_CellChat_candidate_edges.tsv", "S22,S29", "inferred candidate edge", "Focused ANGPT2 CellChat results; computational inference."),
    ("S23_drugreflector_top50.tsv", "13_virtual_perturbation/drugreflector/drugreflector_top50.tsv", "S23", "compound/stratum result", "DrugReflector reverse-signature ranking."),
    ("S25_bulk_validation_programs.tsv", "08_program_validation/GSE62452/GSE62452_paired_program_results.tsv", "S20,S25", "program contrast", "GSE62452 paired program results."),
    ("S25_GSE71729_program_contrasts.tsv", "08_program_validation/GSE71729/GSE71729_program_contrasts.tsv", "S20,S25", "program contrast", "GSE71729 program contrasts."),
    ("S26_GSE297144_anchor.tsv", "09_spatial_validation/GSE297144/GSE297144_patient_spatial_anchor_results.tsv", "S26", "spatial sample", "GSE297144 patient/sample spatial anchor results."),
    ("S27_virtual_ko_top_genes.tsv", "13_virtual_perturbation/closed_loop/supplement_gene_level/EPAS1_KO_top_genes_by_dataset.tsv", "S27", "cohort-gene result", "Virtual EPAS1 knockout gene-level audit."),
    ("S27_virtual_ko_rank_summary.tsv", "13_virtual_perturbation/closed_loop/supplement_gene_level/EPAS1_KO_cross_dataset_rank_summary.tsv", "S27", "gene", "Cross-dataset virtual knockout rank summary."),
    ("S28_independent_localization.tsv", "18_external_strengthening/GSE202051_independent_localization_summary.tsv", "S28", "independent cohort", "GSE202051 endothelial localization summary."),
    ("S30_immune_receiver_statistics.tsv", "20_GSE205049_immune_receiver/GSE205049_immune_receiver_paired_statistics.tsv", "S30", "paired endpoint", "GSE205049 paired immune-receiver sensitivity statistics."),
    ("S30_immune_receiver_patient_data.tsv", "20_GSE205049_immune_receiver/GSE205049_immune_receiver_patient_pseudobulk.tsv", "S30", "patient-endpoint summary", "GSE205049 patient-level immune-receiver values."),
    ("S31_directionality_audit.tsv", "21_revision_strengthening/EPAS1_directionality_audit.tsv", "S31", "evidence row", "Directionality audit across evidence layers."),
    ("S32_conditional_patient_coefficients.tsv", "21_revision_strengthening/EPAS1_spatial_conditional_model_patient_coefficients.tsv", "S32", "patient-outcome coefficient", "Within-patient conditional spatial model coefficients."),
    ("S32_GSE297144_anchor_control.tsv", "21_revision_strengthening/EPAS1_spatial_conditional_model_GSE297144_anchor_control.tsv", "S32", "sample coefficient", "GSE297144 endothelial-anchor control."),
]

SUPPLEMENTARY_TABLES = [
    ("Supplementary_Table_S1_data_resources.tsv", [
        ("public_data_manifest", REPOSITORY_ROOT / "metadata" / "public_data_manifest.tsv"),
    ], "Cohort accession, resource role, and intended inferential unit."),
    ("Supplementary_Table_S2_pseudobulk_and_candidates.tsv", [
        ("pseudobulk_metadata", ANALYSIS_ROOT / "05_pseudobulk/PDAC_patient_pseudobulk_metadata_all.tsv"),
        ("DE_concordance", ANALYSIS_ROOT / "06_pseudobulk_DE/celltype_DE_concordance.tsv"),
        ("target_prioritization", ANALYSIS_ROOT / "08_program_validation/endothelial_target_prioritization.tsv"),
    ], "Patient/sample-level pseudobulk metadata and candidate convergence."),
    ("Supplementary_Table_S3_spatial_and_interface.tsv", [
        ("spatial_anchor", ANALYSIS_ROOT / "09_spatial_validation/GSE282302/GSE282302_patient_spatial_anchor_results.tsv"),
        ("SecAct", ANALYSIS_ROOT / "09_spatial_validation/SecAct/final_evidence_matrix.tsv"),
        ("CellChat", ANALYSIS_ROOT / "17_P0_reinforcement/ANGPT2_CellChat_candidate_edges.tsv"),
        ("conditional_spatial", ANALYSIS_ROOT / "21_revision_strengthening/EPAS1_spatial_conditional_model_summary.tsv"),
    ], "Spatial, SecAct, CellChat, and conditioned-model outputs; inference units are explicit."),
    ("Supplementary_Table_S4_virtual_perturbation.tsv", [
        ("signed_program_responses", ANALYSIS_ROOT / "21_revision_strengthening/EPAS1_virtual_perturbation_signed_program_responses.tsv"),
        ("network_summary", ANALYSIS_ROOT / "13_virtual_perturbation/virtual_perturbation_network_program_summary.tsv"),
    ], "Signed virtual-perturbation proxy responses and network-significance summaries."),
    ("Supplementary_Table_S5_validation_boundaries.tsv", [
        ("MR_instrument_gate", ANALYSIS_ROOT / "11_methodology_review/supplemental_gate/MR_instrument_gate_summary.tsv"),
        ("GSE21501_survival", ANALYSIS_ROOT / "08_program_validation/GSE21501/GSE21501_EPAS1_external_survival.tsv"),
        ("CPTAC_survival", ANALYSIS_ROOT / "11_methodology_review/CPTAC_PAAD/CPTAC_protein_survival_cox.tsv"),
        ("DepMap_PRISM", ANALYSIS_ROOT / "11_methodology_review/DepMap_PRISM_pharmacology/DepMap_CRISPR_PRISM_association_results.tsv"),
        ("survival_ML", ANALYSIS_ROOT / "14_prediction/epas1_survival_ml/survival_model_performance.tsv"),
    ], "MR, survival, protein, dependency/pharmacology, and machine-learning boundary analyses."),
    ("Supplementary_Table_S6_pharmacology_and_docking.tsv", [
        ("DrugReflector", ANALYSIS_ROOT / "13_virtual_perturbation/drugreflector/drugreflector_top50.tsv"),
        ("PRISM", ANALYSIS_ROOT / "11_methodology_review/PRISM20Q2/PRISM20Q2_EPAS1_candidate_drug_results.tsv"),
        ("docking_summary", ANALYSIS_ROOT / "13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_exhaustiveness16_summary.tsv"),
        ("docking_contacts", ANALYSIS_ROOT / "13_virtual_perturbation/docking/EPAS1_multiconformer/final_exhaustiveness16/final_contact_by_conformer.tsv"),
    ], "DrugReflector, PRISM, docking, and geometry-based contact summaries."),
]


def sanitize_table(source: Path, destination: Path) -> tuple[int, int, list[str]]:
    table = pd.read_csv(source, sep="\t", low_memory=False)
    dropped = [column for column in table.columns if PATH_COLUMN.search(str(column))]
    table = table.drop(columns=dropped)
    for column in table.select_dtypes(include="object"):
        table[column] = table[column].map(lambda value: "[local path removed]" if isinstance(value, str) and ABSOLUTE_PATH.search(value) else value)
    table.to_csv(destination, sep="\t", index=False)
    return len(table), len(table.columns), dropped


def checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_combined(output_name: str, parts: list[tuple[str, Path]], description: str) -> dict[str, str | int]:
    frames = []
    dropped_columns = set()
    for label, source in parts:
        if not source.exists():
            raise FileNotFoundError(source)
        frame = pd.read_csv(source, sep="\t", low_memory=False)
        dropped = [column for column in frame.columns if PATH_COLUMN.search(str(column))]
        dropped_columns.update(dropped)
        frame = frame.drop(columns=dropped)
        for column in frame.select_dtypes(include="object"):
            frame[column] = frame[column].map(lambda value: "[local path removed]" if isinstance(value, str) and ABSOLUTE_PATH.search(value) else value)
        frame.insert(0, "source_table", label)
        frames.append(frame)
    destination = OUTPUT_ROOT / output_name
    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined.to_csv(destination, sep="\t", index=False)
    return {
        "file": f"source_data/supplementary/{output_name}",
        "panels": "Supplementary Table",
        "row_unit": "varies by source_table",
        "description": description,
        "upstream_source": "; ".join(
            f"analysis/{source.relative_to(ANALYSIS_ROOT).as_posix()}"
            if source.is_relative_to(ANALYSIS_ROOT)
            else source.relative_to(REPOSITORY_ROOT).as_posix()
            for _, source in parts
        ),
        "rows": len(combined),
        "columns": len(combined.columns),
        "dropped_local_path_columns": ";".join(sorted(dropped_columns)) if dropped_columns else "none",
        "sha256": checksum(destination),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    records = []
    for output_name, relative_source, panels, row_unit, description in EXPORTS:
        source = ANALYSIS_ROOT / relative_source
        if not source.exists():
            raise FileNotFoundError(source)
        destination = OUTPUT_ROOT / output_name
        rows, columns, dropped = sanitize_table(source, destination)
        records.append({
            "file": f"source_data/supplementary/{output_name}",
            "panels": panels,
            "row_unit": row_unit,
            "description": description,
            "upstream_source": f"analysis/{relative_source}",
            "rows": rows,
            "columns": columns,
            "dropped_local_path_columns": ";".join(dropped) if dropped else "none",
            "sha256": checksum(destination),
        })
    for output_name, parts, description in SUPPLEMENTARY_TABLES:
        records.append(write_combined(output_name, parts, description))
    pd.DataFrame(records).to_csv(OUTPUT_ROOT / "manifest.tsv", sep="\t", index=False)
    print(f"Exported {len(records)} supplementary source-data files to {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
