# AGENTS.md

## Project Shape
- This is a small R analysis repo, not a packaged R project: there is no `renv.lock`, `DESCRIPTION`, RStudio project file, CI workflow, or configured test/lint runner.
- Raw/analysis inputs are committed CSVs under `dat/`; current scripts read `rdcp_exhal_analysis.csv`, `dat_ch_enose_20260105.csv`, and `dat_ch_enose_20260808.csv` directly.

## R Workflow
- Run analysis scripts from `analysis/`, because they use relative paths such as `source("./data_prep.R")` and `read.csv("../dat/...")`.
- Main entrypoints are ordered: `setup.R` defines/install-loads packages and helpers, `data_prep.R` sources setup and creates `final_ast_hc_dat`, and `pca_descriptive.R` sources data prep and builds the PCA plot.
- Focused smoke check: `Rscript pca_descriptive.R` from `analysis/`. `Rscript analysis/pca_descriptive.R` from the repo root will not resolve the script's relative paths.
- `setup.R` calls `install.packages()` for missing packages via `packageManage()`. Do not run it casually in constrained/offline environments unless package installation is acceptable.

## Editing Notes
- Preserve the current script style unless making a broader cleanup: base assignment with `=`, section divider comments, and explicit `ggplot2::` qualification only in parts of `setup.R`.
- `data_prep.R` currently has the merged CSV write disabled (`write.csv(..., "../dat/merged_e_dat_rdcp.csv")` is commented); do not create derived data files unless requested.
