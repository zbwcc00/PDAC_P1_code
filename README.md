# EPAS1-associated endothelial remodeling in pancreatic ductal adenocarcinoma

Code and processed, non-identifying result tables supporting the manuscript, *An EPAS1-associated endothelial remodeling signal defines a vessel-centered niche in pancreatic ductal adenocarcinoma*.

## Scope

This repository reproduces the computational workflow used to identify and characterize an EPAS1-associated, vascularly localized transcriptomic signal in pancreatic ductal adenocarcinoma (PDAC). The study integrates public single-cell RNA-sequencing, bulk-expression, spatial-transcriptomic, pharmacology, and structural data. It does not claim that EPAS1 is a causal driver, that a predicted ligand-receptor interaction is direct, or that triclabendazole binds EPAS1 or is therapeutically effective in PDAC.

Raw data are not redistributed here. They remain available from their original public repositories and must be downloaded in accordance with the respective source terms.

## Repository layout

- `analysis/`: analysis scripts retained from the study workflow, organized by stage.
- `config/`: portable configuration templates. Copy the template and set a local project root before running scripts.
- `metadata/`: public-data accession manifest and analysis-to-figure/source-data maps.
- `source_data/`: reserved for small, non-identifying figure-source tables; see `source_data/README.md`.
- `environment/`: R session records and a Python dependency specification.
- `docs/`: workflow, data-access, and reproducibility guidance.
- `scripts/`: repository-level setup and validation helpers.

## Data sources

The main public accessions are GSE154778, GSE155698, GSE212966, GSE62452, GSE71729, GSE282302, GSE297144, GSE327056, GSE202051, GSE205049, and GSE21501. Additional resource-based analyses use TCGA-PAAD, CPTAC-PAAD, DepMap, PRISM, GTEx, the eQTL Catalogue, and the GWAS Catalog. See `metadata/public_data_manifest.tsv` for the role, unit of analysis, and intended use of each source.

## Quick start

1. Create an R environment that includes the packages listed in `environment/R_packages.md` and a Python environment from `environment/requirements.txt`.
2. Download the public datasets listed in `metadata/public_data_manifest.tsv` into a directory outside this repository, preserving their original accession-based names.
3. Set `PDAC_PROJECT_ROOT` to the local study directory, for example in PowerShell: `$env:PDAC_PROJECT_ROOT = 'D:\path\to\PDAC_P1'`.
4. Run the manuscript-critical stages with `python run_pipeline.py --dry-run` to inspect commands, then `python run_pipeline.py` to execute them. The retained non-core scripts still require path migration; see `docs/reproducibility_notes.md`.
5. Compare regenerated summaries with the deposited processed tables and figure source data once the public archive DOI is assigned.

## Source-data verification

To regenerate the compact main-figure and table source data from the study archive:

```powershell
$env:PDAC_PROJECT_ROOT = 'D:\path\to\PDAC_P1'
python scripts/export_source_data.py
python scripts/verify_source_data.py
```

The verification script checks the manuscript's key cohort counts, independent-localization statistics, conditional spatial estimates, pharmacology summary, and removal of local docking paths. Supplementary source-data files require an additional panel-level review before public deposition.

## Reproducibility boundary

Patient/sample is the inferential unit for pseudobulk and spatial-replication analyses. Computational outputs from SecAct, CellChat, virtual perturbation, DrugReflector, PRISM, and docking are hypothesis-generating layers. They should not be interpreted as direct signaling, genetic causality, target engagement, or clinical efficacy.

## Data and code availability

The raw and processed third-party datasets are available from their original repositories under the accession numbers in `metadata/public_data_manifest.tsv`. The versioned analysis code and processed, non-identifying source-data tables are archived in Zenodo: version `v0.1.2`, [10.5281/zenodo.22647158](https://doi.org/10.5281/zenodo.22647158); concept DOI, [10.5281/zenodo.22647157](https://doi.org/10.5281/zenodo.22647157). The matching GitHub release is [v0.1.2](https://github.com/zbwcc00/PDAC_P1_code/releases/tag/v0.1.2).

## Licence

Code is intended for release under the MIT License. Confirm institutional policy and third-party software licences before creating the public release. Do not upload third-party datasets, proprietary software, credentials, or files with personally identifying information.
