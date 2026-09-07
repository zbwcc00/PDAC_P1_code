#!/usr/bin/env python3
"""Export compact, non-identifying source-data tables from the study archive."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pandas as pd

REPOSITORY_ROOT_FROM_SCRIPT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT_FROM_SCRIPT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT_FROM_SCRIPT))

from config.paths import ANALYSIS_ROOT, REPOSITORY_ROOT


SOURCE_DATA_ROOT = REPOSITORY_ROOT / "source_data"
SOURCE_DATA_ROOT.mkdir(parents=True, exist_ok=True)


def copy_tsv(relative_source: str, output_name: str) -> None:
    source = ANALYSIS_ROOT.parent / relative_source
    destination = SOURCE_DATA_ROOT / output_name
    if not source.exists():
        raise FileNotFoundError(source)
    shutil.copyfile(source, destination)


def write_combined(output_name: str, parts: list[tuple[str, str]], source_column: str = "source_table") -> None:
    frames = []
    for label, relative_source in parts:
        source = ANALYSIS_ROOT.parent / relative_source
        frame = pd.read_csv(source, sep="\t")
        frame.insert(0, source_column, label)
        frames.append(frame)
    pd.concat(frames, ignore_index=True, sort=False).to_csv(SOURCE_DATA_ROOT / output_name, sep="\t", index=False)


def export_docking() -> None:
    source = ANALYSIS_ROOT.parent / "analysis/21_revision_strengthening/Figure_5_docking_source_data.tsv"
    frame = pd.read_csv(source, sep="\t")
    frame = frame.drop(columns=["pose_file", "log_file"], errors="ignore")
    frame.insert(0, "source_note", "Absolute local paths removed for public release")
    frame.to_csv(SOURCE_DATA_ROOT / "Figure5_docking_source_data.tsv", sep="\t", index=False)


def main() -> None:
    copy_tsv("analysis/02_scrna_qc/PDAC_scRNA_dataset_summary.tsv", "Figure1_dataset_scale.tsv")
    copy_tsv("analysis/07_candidate_screen/candidate_screen_summary.tsv", "Figure1_candidate_screen_summary.tsv")
    write_combined(
        "Figure2_bulk_spatial_source_data.tsv",
        [
            ("GSE62452_paired_programs", "analysis/08_program_validation/GSE62452/GSE62452_paired_program_results.tsv"),
            ("GSE71729_program_contrasts", "analysis/08_program_validation/GSE71729/GSE71729_program_contrasts.tsv"),
            ("GSE282302_patient_spatial", "analysis/09_spatial_validation/GSE282302/GSE282302_patient_spatial_anchor_results.tsv"),
            ("GSE297144_patient_spatial", "analysis/09_spatial_validation/GSE297144/GSE297144_patient_spatial_anchor_results.tsv"),
        ],
    )
    write_combined(
        "Figure3_signaling_immune_source_data.tsv",
        [
            ("SecAct_cross_cohort", "analysis/09_spatial_validation/SecAct/SecAct_cross_cohort_direction_summary.tsv"),
            ("CellChat_ANGPT2_edges", "analysis/17_P0_reinforcement/ANGPT2_CellChat_candidate_edges.tsv"),
            ("conditional_spatial", "analysis/21_revision_strengthening/EPAS1_spatial_conditional_model_summary.tsv"),
        ],
    )
    copy_tsv("analysis/18_external_strengthening/GSE202051_independent_localization_summary.tsv", "Figure2_GSE202051_localization_summary.tsv")
    copy_tsv("analysis/21_revision_strengthening/EPAS1_virtual_perturbation_signed_program_responses.tsv", "Figure4_virtual_perturbation_source_data.tsv")
    copy_tsv("analysis/13_virtual_perturbation/virtual_perturbation_network_program_summary.tsv", "Figure4_network_program_summary.tsv")
    write_combined(
        "Figure5_pharmacology_source_data.tsv",
        [
            ("PRISM_EPAS1_candidates", "analysis/11_methodology_review/PRISM20Q2/PRISM20Q2_EPAS1_candidate_drug_results.tsv"),
            ("EPAS1_structure_evidence", "analysis/13_virtual_perturbation/docking/EPAS1_multiconformer/epas1_final_structure_evidence.tsv"),
        ],
    )
    export_docking()
    copy_tsv("analysis/02_scrna_qc/PDAC_scRNA_dataset_summary.tsv", "Table1_single_cell_cohort_summary.tsv")
    copy_tsv("analysis/21_revision_strengthening/EPAS1_spatial_conditional_model_summary.tsv", "Table2_conditional_model_results.tsv")
    copy_tsv("analysis/18_external_strengthening/GSE202051_independent_localization_summary.tsv", "Table2_independent_localization.tsv")


if __name__ == "__main__":
    main()
