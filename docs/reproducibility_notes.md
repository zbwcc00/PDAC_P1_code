# Reproducibility notes and pre-release actions

Several retained scripts contain study-era absolute paths such as `D:/PDAC_P1`, `D:/第二篇大论文`, or `T:/`. Before public release, replace these paths with imports from `config/project_paths.R` or `config/project_paths.py`, then perform a clean rerun of the documented main workflow.

The repository currently preserves the original staged scripts for auditability. A full public release should additionally include:

1. a top-level driver (for example, `run_all.R` or Snakemake workflow) that documents dependencies and order;
2. a locked R environment (`renv.lock`) and pinned Python requirements;
3. processed source-data tables mapped to every main-figure panel and every manuscript statistic;
4. a versioned Zenodo release DOI linked to the manuscript; and
5. a repository release tag matching the submitted manuscript version.

Do not upload local RDS/H5AD matrices, credentials, commercial software, or third-party data without confirming redistribution permissions.
