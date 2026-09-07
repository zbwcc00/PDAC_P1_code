#!/usr/bin/env python3
"""Run the ten manuscript-critical analysis stages in dependency order."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from config.paths import ANALYSIS_ROOT, PROJECT_ROOT, REPOSITORY_ROOT


STAGES = [
    ("qc", "R", Path("analysis/02_scrna_qc/run_scrna_qc_seurat.R")),
    ("annotation", "R", Path("analysis/03_scrna_annotation/annotate_scrna_major_lineages.R")),
    ("pseudobulk", "R", Path("analysis/05_pseudobulk/aggregate_patient_pseudobulk.R")),
    ("differential_expression", "R", Path("analysis/06_pseudobulk_DE/run_edgeR_DESeq2_celltype.R")),
    ("spatial_anchor", "R", Path("analysis/09_spatial_validation/run_GSE282302_spatial_anchor.R")),
    ("cellchat", "R", Path("analysis/10_communication_pseudotime/run_cellchat_primary.R")),
    ("virtual_perturbation", "R", Path("analysis/13_virtual_perturbation/02_run_virtual_perturbation.R")),
    ("external_localization", "PYTHON", Path("analysis/18_external_strengthening/01_audit_GSE202051_and_GSE300595.py")),
    ("conditional_spatial", "PYTHON", Path("analysis/21_revision_strengthening/run_directionality_spatial_models_figures.py")),
    ("figure_assembly", "R", Path("analysis/15_priority_figures/38_rebuild_main_figures_unified.R")),
]


def command_for(kind: str, script: Path) -> list[str]:
    executable = "Rscript" if kind == "R" else sys.executable
    return [executable, str(REPOSITORY_ROOT / script)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=[name for name, _, _ in STAGES], action="append")
    parser.add_argument("--from-stage", choices=[name for name, _, _ in STAGES])
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them")
    args = parser.parse_args()

    if not PROJECT_ROOT.exists():
        parser.error(f"PDAC_PROJECT_ROOT does not exist: {PROJECT_ROOT}")
    selected = STAGES
    if args.from_stage:
        start = next(index for index, (name, _, _) in enumerate(STAGES) if name == args.from_stage)
        selected = STAGES[start:]
    if args.stage:
        requested = set(args.stage)
        selected = [stage for stage in STAGES if stage[0] in requested]
    environment = os.environ.copy()
    environment["PDAC_CODE_REPO_ROOT"] = str(REPOSITORY_ROOT)
    environment["PDAC_PROJECT_ROOT"] = str(PROJECT_ROOT)
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, [str(REPOSITORY_ROOT), environment.get("PYTHONPATH", "")]))

    print(f"Repository: {REPOSITORY_ROOT}")
    print(f"Project data root: {PROJECT_ROOT}")
    print(f"Analysis root: {ANALYSIS_ROOT}")
    for name, kind, script in selected:
        command = command_for(kind, script)
        print(f"\n[{name}] {' '.join(command)}")
        if args.dry_run:
            continue
        result = subprocess.run(command, cwd=REPOSITORY_ROOT, env=environment)
        if result.returncode:
            print(f"Stage failed: {name} (exit {result.returncode})", file=sys.stderr)
            return result.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
