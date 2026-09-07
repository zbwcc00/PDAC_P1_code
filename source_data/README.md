# Source data

This directory contains compact, non-identifying TSV files for the main figures and tables. They were exported on 2026-09-07 with `scripts/export_source_data.py`. Third-party raw matrices are not redistributed here; see `metadata/public_data_manifest.tsv`.

Current release files include Figure 1 cohort and candidate summaries; Figure 2 bulk, spatial, and GSE202051 localization tables; Figure 3 signaling and conditional-model outputs; Figure 4 virtual-perturbation and network summaries; Figure 5 pharmacology and docking tables; and Table 1-2 summaries.

Each table should provide:

- a stable column name for each variable;
- the biological unit and accession represented by each row;
- units and missing-value codes;
- the analysis script that generated the value; and
- a short panel-level description. The row unit is indicated by the source table and column names (cohort, patient/sample, spatial endpoint, program, or docking structure).

## Conventions

- Files are UTF-8, tab-delimited TSV; missing values are empty and are not zeroes.
- `P.Value`, `adj.P.Val`, `*_p`, and `*_q` are two-sided unless the manuscript specifies otherwise.
- `rho` denotes Spearman correlation; `logFC`/`lfc` denote method-defined log fold change; affinity values are kcal/mol.
- Local docking pose and log paths were removed from the public table. SecAct, CellChat, virtual perturbation, and docking outputs remain computational/inferred evidence as described in the manuscript.

After changing the manuscript or analysis, rerun the export script, update `../metadata/figure_source_data_map.tsv`, run the repository audit, and verify that the tables reproduce the submitted manuscript.
