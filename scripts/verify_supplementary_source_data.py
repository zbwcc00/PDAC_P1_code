#!/usr/bin/env python3
"""Verify the integrity and panel coverage of supplementary source data."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENTARY = ROOT / "source_data" / "supplementary"
MANIFEST = SUPPLEMENTARY / "manifest.tsv"
ABSOLUTE_PATH = re.compile(r"(?:(?:[A-Za-z]:\\)|(?:/Users/)|(?:/home/))")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    manifest = pd.read_csv(MANIFEST, sep="\t")
    expected_figures = {f"S{index}" for index in range(1, 33)}
    covered_figures = set()
    failures: list[str] = []

    for record in manifest.itertuples(index=False):
        path = ROOT / record.file
        if not path.exists():
            failures.append(f"missing file: {record.file}")
            continue
        if sha256(path) != record.sha256:
            failures.append(f"checksum mismatch: {record.file}")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        values = frame.select_dtypes(include="object").fillna("").astype(str)
        if values.map(lambda value: bool(ABSOLUTE_PATH.search(value))).any().any():
            failures.append(f"absolute local path found: {record.file}")
        covered_figures.update(re.findall(r"S(\d+)", str(record.panels)))

    missing_figures = expected_figures - {f"S{number}" for number in covered_figures}
    if missing_figures:
        failures.append("missing figure coverage: " + ", ".join(sorted(missing_figures, key=lambda value: int(value[1:]))))
    expected_tables = {f"Supplementary_Table_S{index}_" for index in range(1, 7)}
    names = {Path(value).name for value in manifest.file}
    for prefix in expected_tables:
        if not any(name.startswith(prefix) for name in names):
            failures.append(f"missing supplementary table: {prefix[:-1]}")

    if failures:
        for failure in failures:
            print(f"FAIL\t{failure}")
        return 1
    print(f"PASS\t{len(manifest)} source-data tables present with matching checksums")
    print("PASS\tSupplementary Figures S1-S32 covered")
    print("PASS\tSupplementary Tables S1-S6 covered")
    print("PASS\tNo absolute local paths detected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
