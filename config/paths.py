"""Portable paths for the PDAC P1 analysis repository."""

from __future__ import annotations

import os
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(os.environ.get("PDAC_PROJECT_ROOT", REPOSITORY_ROOT)).expanduser().resolve()
RAW_DATA_ROOT = PROJECT_ROOT / "data"
EXTERNAL_DATA_ROOT = PROJECT_ROOT / "data_external"
ANALYSIS_ROOT = PROJECT_ROOT / "analysis"
RESULTS_ROOT = PROJECT_ROOT / "results"
FIGURE_ROOT = ANALYSIS_ROOT / "15_priority_figures"


def describe_paths() -> dict[str, Path]:
    return {
        "repository_root": REPOSITORY_ROOT,
        "project_root": PROJECT_ROOT,
        "raw_data_root": RAW_DATA_ROOT,
        "external_data_root": EXTERNAL_DATA_ROOT,
        "analysis_root": ANALYSIS_ROOT,
        "results_root": RESULTS_ROOT,
        "figure_root": FIGURE_ROOT,
    }
