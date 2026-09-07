# Workflow

## Primary analysis stages

The ten manuscript-critical stages can be listed without running them with:

```powershell
python run_pipeline.py --dry-run
```

After `PDAC_PROJECT_ROOT` is set and dependencies/data are available, run all stages with `python run_pipeline.py`, or resume from a named stage with `python run_pipeline.py --from-stage spatial_anchor`.

1. `analysis/02_scrna_qc/`: parse GEO manifests, perform sample-aware quality control, and retain QC reports.
2. `analysis/03_scrna_annotation/`: normalize retained objects and assign conservative major-lineage labels.
3. `analysis/04_cnv_copykat/`: generate exploratory CopyKAT/CNV support for malignant epithelial annotations.
4. `analysis/05_pseudobulk/` and `analysis/06_pseudobulk_DE/`: aggregate counts at patient/sample level and perform lineage-restricted edgeR/DESeq2 comparisons.
5. `analysis/08_program_validation/` through `analysis/12_functional_immune/`: validate program behavior in bulk and spatial cohorts and run activity/communication inference.
6. `analysis/13_virtual_perturbation/` and `analysis/14_prediction/`: perform virtual network perturbation and boundary/sensitivity analyses.
7. `analysis/18_external_strengthening/`, `analysis/20_GSE205049_immune_receiver/`, and `analysis/21_revision_strengthening/`: run independent localization, receiver-side, and conditioned spatial analyses.
8. `analysis/15_priority_figures/`: assemble figures from finalized summaries.

## Important restrictions

- GSE154778 primary-versus-metastasis comparisons are exploratory and are not normal-tissue controls.
- GSE155698 adjacent-normal endothelial comparisons have limited replication and must not be used alone as core evidence.
- CopyKAT/CNV calls are exploratory support, not definitive normal/malignant labels.
- The VM scripts in `analysis/22_vm_exploratory/` are audit/exploratory work and are not part of the current manuscript's confirmatory claim set.
- Do not pool spatial spots as independent biological replicates; calculate within-patient/sample estimates first.
