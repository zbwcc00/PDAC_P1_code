# Supplementary source data

This directory contains 45 compact, non-identifying TSV result tables that support Supplementary Figures S1-S32 and Supplementary Tables S1-S6 of the manuscript. It was generated with `scripts/export_supplementary_source_data.py` on 2026-09-07.

`manifest.tsv` is the authoritative map: it gives each file's manuscript panel(s), row unit, description, upstream analysis result, dimensions, removed local-path columns, and SHA-256 checksum. Use the manifest rather than inferring the intended analysis from a filename.

## Scope and conventions

- These are processed result tables, not redistribution copies of GEO, TCGA, CPTAC, DepMap, PRISM, or other third-party raw datasets. Original accessions and resource roles are listed in `../../metadata/public_data_manifest.tsv`.
- Empty fields mean that a variable does not apply to a source-table row, rather than zero.
- Statistical outputs retain the definitions established by their upstream analyses. Patient/sample is the inferential unit where stated in the manuscript; spot- and cell-level data are used only for descriptive displays or computational summaries.
- SecAct, CellChat, virtual perturbation, DrugReflector, and docking tables are computational/inferred outputs. They do not establish direct signaling, causality, target engagement, or efficacy.
- Columns that name a local file, log, or directory are removed. Any remaining value resembling an absolute local path is replaced with `[local path removed]`.

## Supplementary tables

`Supplementary_Table_S1_data_resources.tsv` through `Supplementary_Table_S6_pharmacology_and_docking.tsv` correspond exactly to manuscript Supplementary Tables S1-S6. Their `source_table` column identifies the contributing analysis table; other blank columns are not applicable to that table type.

## Regeneration

Set `PDAC_PROJECT_ROOT` to the local study directory that contains the original analysis outputs, then run:

```powershell
python scripts/export_supplementary_source_data.py
python scripts/verify_supplementary_source_data.py
```

Re-run both commands after modifying an upstream analysis or figure. The checksums in `manifest.tsv` are regenerated with the tables and should not be edited manually.
