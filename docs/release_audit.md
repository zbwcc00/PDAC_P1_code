# Public-release audit

## Current status

The repository is a **working draft**, not yet a public release. The ten manuscript-critical scripts are portable and the compact main-figure/table source-data package has been exported and checked. The wider historical script archive still contains study-era absolute paths, so it is not yet a one-command reproducible release.

## Release gate

Run the following from the repository root before public upload:

```powershell
python scripts/audit_repository.py
```

The audit must report zero absolute paths in the ten core scripts and zero unreviewed credential-pattern hits. `python scripts/verify_source_data.py` must pass before marking main source-data files `READY_FOR_RELEASE`. Before calling the full historical script archive one-command reproducible, migrate the remaining scripts too. Review all generated figures and confirm that the exact released scripts regenerate the deposited source-data tables.

## Priority migration set

Convert these scripts first, as they support the manuscript's central claims:

1. `analysis/02_scrna_qc/run_scrna_qc_seurat.R`
2. `analysis/03_scrna_annotation/annotate_scrna_major_lineages.R`
3. `analysis/05_pseudobulk/aggregate_patient_pseudobulk.R`
4. `analysis/06_pseudobulk_DE/run_edgeR_DESeq2_celltype.R`
5. `analysis/09_spatial_validation/run_GSE282302_spatial_anchor.R`
6. `analysis/10_communication_pseudotime/run_cellchat_primary.R`
7. `analysis/13_virtual_perturbation/02_run_virtual_perturbation.R`
8. `analysis/18_external_strengthening/01_audit_GSE202051_and_GSE300595.py`
9. `analysis/21_revision_strengthening/run_directionality_spatial_models_figures.py`
10. `analysis/15_priority_figures/38_rebuild_main_figures_unified.R`

Each migrated script should import `config/paths.R` or `config/paths.py`, resolve all inputs from the configured project root, and document its input files and output files. Existing scripts retain their analysis-stage output directories so that their results remain compatible with the manuscript audit trail; a later release may redirect generated outputs to `RESULTS_ROOT`.

## Source-data deposit

The `source_data/` directory now contains small, non-identifying TSV files for the main figures and tables, with a panel map in `metadata/figure_source_data_map.tsv`. Run `python scripts/verify_source_data.py` after every export. Supplementary source-data packages remain `TO_DEPOSIT` until their panels are mapped and checked. Raw GEO matrices should not be duplicated in this repository.

## Zenodo release

The GitHub release `v0.1.0` is available at `https://github.com/zbwcc00/PDAC_P1_code/releases/tag/v0.1.0`. Authorize GitHub in Zenodo, enable this repository, and use that release to create the archival record. Verify the draft record's creators, affiliation, description, keywords, version, related GitHub URL, licence, file manifest, and public visibility before publishing. Add the final DOI to the manuscript's Code Availability statement and Data Availability statement if source-data tables are archived in the same record. Do not mint or report a DOI until Zenodo assigns it.
