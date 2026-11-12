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

## Interpretation
- Lambda selection and model fitting use training data only. The fitted model is then evaluated on the held-out test batch.
- The reported AUC and other performance measures describe predictions made using all 13 sensors together.

### Why coefficient size does not directly rank sensor importance

`standardize = TRUE` standardizes the sensor values internally during model fitting. However, `coef(fit_final, s = lambda_best)` converts the fitted coefficients back to the original sensor scales. This is normal glmnet behaviour.

Each returned coefficient describes the change in the model's asthma score (log-odds) for a one-unit increase in that sensor, holding the other sensor values fixed. It is not a percentage-point change in predicted probability.

Sensors vary over different numerical ranges. A sensor with a narrow range can therefore have a large coefficient: the effect on the model's score depends on **the coefficient multiplied by the change in sensor value**.

For example, S6BH's coefficient of approximately 19.15 does not automatically make it more important than S7, whose coefficient is approximately 10.27. The amount each sensor varies also matters. The raw coefficient sizes alone cannot establish which sensor contributes more to test-set discrimination.

## Optional sensor interpretation
Supplementary permutation-importance and Shapley analyses are available in `analysis_ast_after_filter/ridge_train_test_with_shap_perm_etc.R`. They explore how the fitted model uses individual sensors and complement the main train/test performance analysis.
