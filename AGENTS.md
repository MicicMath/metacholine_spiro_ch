# AGENTS.md

## Project Shape
- This is a small R analysis repo, not a packaged R project: there is no `renv.lock`, `DESCRIPTION`, RStudio project file, CI workflow, or configured test/lint runner.
- Raw/analysis inputs are committed CSVs under `dat/`; current scripts read `rdcp_exhal_analysis.csv`, `dat_ch_enose_20260105.csv`, and `dat_ch_enose_20260808.csv` directly.
- `rdcp_exhal_analysis.csv` is the parsed REDCap export. It is merged with eNose breath sensor readings by measurement ID, not patient ID: REDCap `measurement_id_exhal` is renamed to `m_id`, and the second eNose file's `id` column is renamed to `m_id`.
- Patients can have multiple visits/measurements. Treat `m_id` as the measurement-level key and `patient_id` as repeated subject metadata.

## R Workflow
- Run analysis scripts from `analysis/`, because they use relative paths such as `source("./data_prep.R")` and `read.csv("../dat/...")`.
- Main setup flow: `setup.R` defines/install-loads packages and helpers, and `data_prep.R` sources setup, defines `sensors`/`cols`, imports REDCap plus eNose CSVs, calculates `fev1_z`/`fvc_z`, subsets to HC/asthma/CF diagnoses, and creates `final_ast_hc_cf_dat`.
- `final_ast_hc_cf_dat` is first-visit only, keeps HC plus confirmed disease (`diagnosis_status == 1`), removes non-allergic asthma, and adds `diagnosis_simple` plus `metacholine_response`.
- Analysis entrypoints currently include `pca_descriptive.R`, `univariate_comparison.R`, `incl_cf_univariate_comparison.R`, `correlation_plot.R`, `lasso.R`, `svm_rf.R`, and `elanet.R`; all source `data_prep.R` and start from `final_ast_hc_cf_dat`.
- `pca_descriptive.R`, `univariate_comparison.R`, and `lasso.R` drop CF downstream and compare allergic asthma with Mid/High metacholine response against HC. `incl_cf_univariate_comparison.R` keeps CF and compares AST/CF/HC after the same asthma metacholine filtering.
- `correlation_plot.R` shows pairwise correlations among sensor readings after dropping CF and keeping HC plus asthma with Mid/High metacholine response.
- `lasso.R` uses weighted LASSO logistic regression (`glmnet`, `alpha = 1`) with outer repeated stratified 5-fold CV, inner stratified 5-fold `foldid` for `lambda.1se` selection, held-out fold metrics, feature selection frequency across outer-loop fits, and a vertically averaged ROC plot. Optional permutation testing, ridge sensitivity analysis, and plot export blocks are present but commented by default.
- `svm_rf.R` runs SVM and random forest classifiers to check robustness of the AST-vs-HC eNose signal across classifier families.
- `elanet.R` is an elastic-net extension/sensitivity analysis for AST vs HC classification, tuning `alpha` over `0` to `1` with the same repeated stratified 5-fold outer split structure.
- Focused smoke checks from `analysis/`: `Rscript pca_descriptive.R`, `Rscript univariate_comparison.R`, `Rscript incl_cf_univariate_comparison.R`, `Rscript correlation_plot.R`, `Rscript lasso.R`, `Rscript svm_rf.R`, and `Rscript elanet.R`. Running these as `Rscript analysis/<script>.R` from the repo root will not resolve the scripts' relative paths.
- `setup.R` calls `install.packages()` for missing packages via `packageManage()`, including `glmnet`, `ranger`, and `e1071` for classifier scripts. Do not run it casually in constrained/offline environments unless package installation is acceptable.

## Editing Notes
- Preserve the current script style unless making a broader cleanup: base assignment with `=`, section divider comments, and explicit `ggplot2::` qualification only in parts of `setup.R`.
- Do not create derived data files or plot exports unless requested; plotting save calls should remain commented by default.
- Avoid changing the REDCap/eNose merge from measurement-level `m_id` to patient-level `patient_id` unless explicitly requested.
