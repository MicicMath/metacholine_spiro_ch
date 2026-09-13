# AGENTS.md

## Project Shape
- This is a small R analysis repo, not a packaged R project: there is no `renv.lock`, `DESCRIPTION`, RStudio project file, CI workflow, or configured test/lint runner.
- Raw/analysis inputs are committed CSVs under `dat/`; current scripts read `rdcp_exhal_analysis.csv`, `dat_ch_enose_20260105.csv`, and `dat_ch_enose_20260808.csv` directly.
- `rdcp_exhal_analysis.csv` is the parsed REDCap export. It is merged with eNose breath sensor readings by measurement ID, not patient ID: REDCap `measurement_id_exhal` is renamed to `m_id`, and the second eNose file's `id` column is renamed to `m_id`.
- Patients can have multiple visits/measurements. Treat `m_id` as the measurement-level key and `patient_id` as repeated subject metadata.
- The current main analysis folder is `analysis_ast_after_filter/`. The older `analysis/` folder remains in the repo as legacy/contextual work and should not be assumed to be the active workflow unless explicitly requested.

## R Workflow
- Run current analysis scripts from `analysis_ast_after_filter/`, because they use relative paths such as `source("./data_prep.R")` and `read.csv("../dat/...")`. Running them as `Rscript analysis_ast_after_filter/<script>.R` from the repo root will not resolve the scripts' relative paths.
- Current setup flow: `setup.R` defines/install-loads packages and helpers; `source_selected_ast_hc.R` defines manually selected asthma and healthy-control measurement IDs; `data_prep.R` sources both files, imports REDCap plus eNose CSVs, filters REDCap to the selected measurement IDs, adds eNose batch labels, calculates `fev1_z`/`fvc_z`, subsets to HC plus asthma diagnoses, and creates `final_ast_hc_cf_dat`.
- `final_ast_hc_cf_dat` is the active cross-sectional AST-vs-HC analysis object. The name is historical: current `analysis_ast_after_filter/data_prep.R` comments out CF/PCD, includes allergic and non-allergic asthma, labels both as `diagnosis_simple == "AST"`, keeps healthy controls as `HC`, and filters to first visit (`visit_exhal == 1`).
- The current diagnosis-status filter keeps rows where `diagnosis_status` is `0` or `1`, plus HC rows; a stricter `diagnosis_status == 1` line exists but is commented. Do not change this filter without an explicit analysis decision.
- `metacholine_response` is still derived from `metacholine_test_result`, but the current `analysis_ast_after_filter/` entrypoints generally compare selected AST vs HC directly rather than applying the older Mid/High metacholine-response filtering.
- Current entrypoints include `simple_pca.R`, `univariate_comparison.R`, `correlation_plot.R`, `traces.R`, `lasso.R`, `ridge.R`, `elastic_net.R`, `pca_ridge.R`, `permutation_test.R`, `permutation_test_ridge.R`, and `svm_rf_lasso_dry.R`; all source `data_prep.R` and start from `final_ast_hc_cf_dat`.
- `univariate_comparison.R` inspects HC batch effects, runs AST-vs-HC sensor-level Wilcoxon and Welch t-tests with multiple-testing adjustment/effect sizes, and builds boxplots/radar/spider plots. Table and plot exports are commented by default.
- `correlation_plot.R` builds sensor correlation plots for all samples and by AST/HC group, then runs a permutation test comparing AST and HC correlation matrices (`B = 10000` by default).
- `traces.R` plots sensor intensity over measurement date and fits per-sensor linear models adjusted for diagnosis, time, and operator.
- `lasso.R`, `ridge.R`, `elastic_net.R`, and `pca_ridge.R` run weighted logistic models on AST vs HC using repeated stratified 5-fold outer CV (`R = 100` by default), inner CV for tuning, held-out fold metrics, and ROC summaries. `elastic_net.R` tunes `alpha` over `0, 0.25, 0.5, 0.75, 1`; `pca_ridge.R` applies PCA preprocessing inside the caret training workflow.
- `svm_rf_lasso_dry.R` runs LASSO, SVM, and random forest classifiers in one script to check robustness of the AST-vs-HC eNose signal across classifier families and combines formatted performance summaries in `results_all`.
- `permutation_test.R` and `permutation_test_ridge.R` are heavier permutation workflows (`1000` permutations by default). Avoid running them casually unless runtime and CPU use are acceptable.
- Focused smoke checks from `analysis_ast_after_filter/`: `Rscript simple_pca.R`, `Rscript univariate_comparison.R`, `Rscript correlation_plot.R`, `Rscript traces.R`, `Rscript lasso.R`, `Rscript ridge.R`, `Rscript elastic_net.R`, `Rscript pca_ridge.R`, and `Rscript svm_rf_lasso_dry.R`. Run permutation scripts separately when specifically needed.
- `setup.R` calls `install.packages()` for missing packages via `packageManage()`, including modeling/plotting packages such as `glmnet`, `ranger`, `e1071`, `caret`, `Boruta`, `ggpubr`, and `ROCR`. Do not run it casually in constrained/offline environments unless package installation is acceptable.

## Editing Notes
- Preserve the current script style unless making a broader cleanup: base assignment with `=`, section divider comments, and explicit `ggplot2::` qualification only in parts of `setup.R`.
- Do not create derived data files or plot exports unless requested; plotting save calls should remain commented by default.
- Avoid changing the REDCap/eNose merge from measurement-level `m_id` to patient-level `patient_id` unless explicitly requested.
- Avoid changing the selected AST/HC measurement IDs in `analysis_ast_after_filter/source_selected_ast_hc.R` unless explicitly requested.
