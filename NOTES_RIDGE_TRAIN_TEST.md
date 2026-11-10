# Notes: ridge train/test (2026-11-10)

Script: `analysis_ast_after_filter/ridge_train_test.R` (run from inside `analysis_ast_after_filter/`).

## Data split (from `data_prep.R`)
- **Train**: `final_ast_hc_cf_dat` = batches 20260105 + 20260808, n = 95 (42 AST, 53 HC).
- **Test**: `test_set_final_ast_hc_cf_dat` = batch 20261027, n = 69 (36 AST, 33 HC).
- No patients shared between train and test. Operators JST, SVO and VSH appear only in test.
- Features: the 13 sensors (`sensors`), raw values, no batch adjustment.

## Tuning (training data only)
- Ridge logistic regression (`alpha = 0`) with class-balanced weights (as in `ridge.R`).
- **Repeated 10-fold CV:** 100 × `cv.glmnet`, each repeat with a new stratified 10-fold split.
- **Measure:** binomial deviance. AUC is not possible: n = 95 with 10 folds gives < 10 obs per fold, so glmnet switches to deviance.
- **Lambda grid:** chosen by glmnet. It is identical in every repeat because it is built from the full training data. We set `lambda.min.ratio = 1e-6`, because the default grid (1e-4) ends at 0.040 and the optimum is below that.
- **Selection:** average the 100 CV deviance curves, then take the lambda with the lowest mean. This is lambda.min logic.
  - lambda.1se is not possible with repeated CV, because it needs the SE from a single CV run.
  - Source: `?cv.glmnet` Details: "running cv.glmnet many times, and averaging the error curves". Same as caret `repeatedcv` with `selectionFunction = "best"`.
- **Result:** lambda = 0.0113, inside the grid.

## Final model and test
- Refit `glmnet` on all of the training data (full lambda path), read off at lambda = 0.0113.
- Predict the test set **once**, threshold 0.5.
- Results, with 95% stratified bootstrap CIs (B = 2000):

| Metric | Value (95% CI) |
|---|---|
| AUC | 0.93 (0.85, 1.00) |
| Balanced accuracy | 92.8% (86.7, 98.5) |
| Sensitivity | 91.7% (83.3, 100.0) |
| Specificity | 93.9% (84.8, 100.0) |

- Confusion matrix: 33/36 AST and 31/33 HC correct, so 5 misclassified.

## Caveats
- The test set is now used. If the model or tuning changes because of this result, the test is no longer independent.
- Coefficients are on each sensor's original scale, so their sizes can't be compared between sensors.
- Training batch 1 is mostly HC (8 AST, 29 HC), so batch and diagnosis are partly confounded in training.

## Next round (deferred)
- Permutation feature importance on the test set: shuffle one sensor in `X_test` (1000×) and measure the AUC drop, with an interval.
- Grouped permutation for the correlated sensors S1BH, S3BH and S4BH (Spearman 0.76 to 0.84 in train).
- SHAP: for a linear model it is exact, beta_j × (x_j − training mean). Mainly useful to explain the 5 misclassified test samples.
