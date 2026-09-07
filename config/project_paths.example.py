"""Copy to project_paths.py, edit paths, and keep the local copy out of Git."""
from pathlib import Path

PDAC_PROJECT_ROOT = Path(r"D:/path/to/PDAC_P1")
RAW_DATA_ROOT = PDAC_PROJECT_ROOT / "data"
EXTERNAL_DATA_ROOT = PDAC_PROJECT_ROOT / "data_external"
ANALYSIS_ROOT = PDAC_PROJECT_ROOT / "analysis"
RESULTS_ROOT = PDAC_PROJECT_ROOT / "results"
