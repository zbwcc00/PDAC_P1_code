#!/usr/bin/env python3
"""Audit the public-release repository without reading raw data files."""

from __future__ import annotations

import re
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_SUFFIXES = {".R", ".py", ".Rmd"}
ABSOLUTE_PATH = re.compile(r"(?:[A-Za-z]:[\\/]|/home/|/Users/)")
SENSITIVE_PATTERN = re.compile(
    r"(?:password|api[_ -]?key|authorization\s*:|private[_ -]?key)", re.IGNORECASE
)
CORE_SCRIPTS = [
    Path("analysis/02_scrna_qc/run_scrna_qc_seurat.R"),
    Path("analysis/03_scrna_annotation/annotate_scrna_major_lineages.R"),
    Path("analysis/05_pseudobulk/aggregate_patient_pseudobulk.R"),
    Path("analysis/06_pseudobulk_DE/run_edgeR_DESeq2_celltype.R"),
    Path("analysis/09_spatial_validation/run_GSE282302_spatial_anchor.R"),
    Path("analysis/10_communication_pseudotime/run_cellchat_primary.R"),
    Path("analysis/13_virtual_perturbation/02_run_virtual_perturbation.R"),
    Path("analysis/18_external_strengthening/01_audit_GSE202051_and_GSE300595.py"),
    Path("analysis/21_revision_strengthening/run_directionality_spatial_models_figures.py"),
    Path("analysis/15_priority_figures/38_rebuild_main_figures_unified.R"),
]


def audit() -> int:
    scripts = sorted(
        path for path in (REPOSITORY_ROOT / "analysis").rglob("*")
        if path.is_file() and path.suffix in SCRIPT_SUFFIXES
    )
    absolute_path_scripts: list[Path] = []
    sensitive_hits: list[tuple[Path, int]] = []

    for script in scripts:
        text = script.read_text(encoding="utf-8", errors="replace")
        if ABSOLUTE_PATH.search(text):
            absolute_path_scripts.append(script)
        for line_number, line in enumerate(text.splitlines(), start=1):
            if SENSITIVE_PATTERN.search(line):
                sensitive_hits.append((script, line_number))

    core_absolute_paths = [path for path in CORE_SCRIPTS if ABSOLUTE_PATH.search((REPOSITORY_ROOT / path).read_text(encoding="utf-8", errors="replace"))]
    print(f"Scripts scanned: {len(scripts)}")
    print(f"Scripts with absolute paths (all retained scripts): {len(absolute_path_scripts)}")
    for script in absolute_path_scripts:
        print(f"  PATH  {script.relative_to(REPOSITORY_ROOT)}")
    print(f"Core scripts with absolute paths: {len(core_absolute_paths)}")
    for script in core_absolute_paths:
        print(f"  CORE-PATH {script}")
    print(f"Potential credential-pattern hits: {len(sensitive_hits)}")
    for script, line_number in sensitive_hits:
        print(f"  CHECK {script.relative_to(REPOSITORY_ROOT)}:{line_number}")

    return 1 if core_absolute_paths or sensitive_hits else 0


if __name__ == "__main__":
    sys.exit(audit())
