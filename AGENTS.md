# AGENTS.md

## Project Shape
- This is a small R analysis repo, not a packaged R project: there is no `renv.lock`, `DESCRIPTION`, RStudio project file, CI workflow, or configured test/lint runner.
- Raw/analysis inputs are committed CSVs under `dat/`; current scripts read `rdcp_exhal_analysis.csv`, `dat_ch_enose_20260105.csv`, and `dat_ch_enose_20260808.csv` directly.
- `rdcp_exhal_analysis.csv` is the parsed REDCap export. It is merged with eNose breath sensor readings by measurement ID, not patient ID: REDCap `measurement_id_exhal` is renamed to `m_id`, and the second eNose file's `id` column is renamed to `m_id`.
- Patients can have multiple visits/measurements. Treat `m_id` as the measurement-level key and `patient_id` as repeated subject metadata.

## R Workflow
- Run analysis scripts from `analysis/`, because they use relative paths such as `source("./data_prep.R")` and `read.csv("../dat/...")`.
- Main setup flow: `setup.R` defines/install-loads packages and helpers, and `data_prep.R` sources setup, defines `sensors`/`cols`, imports REDCap plus eNose CSVs, subsets to HC/asthma/CF diagnoses, and creates `final_ast_hc_cf_dat`.
- `final_ast_hc_cf_dat` is first-visit only, keeps HC plus confirmed disease (`diagnosis_status == 1`), removes non-allergic asthma, and adds `diagnosis_simple` plus `metacholine_response`.
- Analysis entrypoints currently include `pca_descriptive.R`, `univariate_comparison.R`, and `incl_cf_univariate_comparison.R`; all source `data_prep.R` and start from `final_ast_hc_cf_dat`.
- `pca_descriptive.R` and `univariate_comparison.R` drop CF downstream and compare allergic asthma with Mid/High metacholine response against HC. `incl_cf_univariate_comparison.R` keeps CF and compares AST/CF/HC after the same asthma metacholine filtering.
- Focused smoke checks from `analysis/`: `Rscript pca_descriptive.R`, `Rscript univariate_comparison.R`, and `Rscript incl_cf_univariate_comparison.R`. Running these as `Rscript analysis/<script>.R` from the repo root will not resolve the scripts' relative paths.
- `setup.R` calls `install.packages()` for missing packages via `packageManage()`. Do not run it casually in constrained/offline environments unless package installation is acceptable.

## Editing Notes
- Preserve the current script style unless making a broader cleanup: base assignment with `=`, section divider comments, and explicit `ggplot2::` qualification only in parts of `setup.R`.
- `data_prep.R` currently has the merged CSV write disabled (`write.csv(..., "../dat/merged_e_dat_rdcp.csv")` is commented); do not create derived data files unless requested.
- Avoid changing the REDCap/eNose merge from measurement-level `m_id` to patient-level `patient_id` unless explicitly requested.
