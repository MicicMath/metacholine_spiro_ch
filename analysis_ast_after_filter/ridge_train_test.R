source("./data_prep.R")

auc_fn = function(y_true, y_prob) {
  ROCR::performance(
    ROCR::prediction(as.numeric(y_prob), y_true),
    "auc"
  )@y.values[[1]]
}

perm_importance_fn = function(model, X_test, y_test, sensors, baseline_auc, lambda_value, B = 1000) {
  perm_list = parallel::mclapply(sensors, function(sensor) {
    drops = replicate(B, {
      X_perm = X_test
      X_perm[, sensor] = sample(X_perm[, sensor])

      y_prob_perm = predict(
        model,
        newx = X_perm,
        s = lambda_value,
        type = "response"
      )

      auc_perm = auc_fn(y_true = y_test, y_prob = y_prob_perm)
      auc_drop = baseline_auc - auc_perm

      return(auc_drop)
    })

    data.frame(
      sensor = sensor,
      mean_auc_drop = mean(drops, na.rm = TRUE),
      sd_auc_drop = sd(drops, na.rm = TRUE),
      ci_low = quantile(drops, probs = 0.025, na.rm = TRUE),
      ci_high = quantile(drops, probs = 0.975, na.rm = TRUE),
      row.names = NULL
    )
  }, mc.cores = max(1, parallel::detectCores() - 2))

  perm_imp = do.call(rbind, perm_list)
  perm_imp = perm_imp[order(-perm_imp$mean_auc_drop), ]

  return(perm_imp)
}

perm_importance_plot_fn = function(perm_imp, xlim) {
  plot_dat = perm_imp
  plot_dat$sensor = factor(
    plot_dat$sensor,
    levels = plot_dat$sensor[order(plot_dat$mean_auc_drop)]
  )

  p = ggplot(
    plot_dat,
    aes(x = mean_auc_drop, y = sensor)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      width = 0,
      orientation = "y"
    ) +
    geom_point(size = 2.5, color = "darkblue") +
    coord_cartesian(xlim = xlim) +
    labs(
      x = "Mean AUC drop after permutation",
      y = "Sensor"
    ) +
    theme_classic()

  return(p)
}

#############################################
# train (batches 1 + 2) and test (batch 3) data
#############################################

train_dat = final_ast_hc_cf_dat
test_dat = test_set_final_ast_hc_cf_dat

table(train_dat$diagnosis_simple)
table(test_dat$diagnosis_simple)

X_train = as.matrix(train_dat[, sensors])
y_train = ifelse(train_dat$diagnosis_simple == "AST", 1, 0)

X_test = as.matrix(test_dat[, sensors])
y_test = ifelse(test_dat$diagnosis_simple == "AST", 1, 0)

###############################
# class-balanced weights (train)
###############################

r = 0.5 # 0.5 is standard balance
w_train = ifelse(
  y_train == 1,
  r * (1 / sum(y_train == 1)),
  (1 - r) * (1 / sum(y_train == 0))
)
w_train = w_train * (length(y_train) / sum(w_train))

##########################
# stratified 10-fold split
##########################

make_foldid = function(y, k = 10) {
  folds = caret::createFolds(
    factor(y, levels = c(1, 0)),
    k = k,
    returnTrain = FALSE
  )

  foldid = rep(0, length(y))

  for (i in seq_along(folds)) {
    foldid[folds[[i]]] = i
  }

  return(foldid)
}

######################################
# repeated 10-fold CV (R x cv.glmnet)
######################################

# 10 folds with n = 95 leaves < 10 obs per fold, so AUC is not possible in
# cv.glmnet; binomial deviance is used. lambda.min.ratio = 1e-6 extends the
# glmnet lambda grid downwards, the default (1e-4) ends above the optimum.
# Averaging the CV error curves over repeats follows ?cv.glmnet (Details).

R = 100

set.seed(101)
cv_runs = lapply(
  1:R,
  function(i) {
    glmnet::cv.glmnet(
      x = X_train,
      y = y_train,
      family = "binomial",
      alpha = 0,
      weights = w_train,
      foldid = make_foldid(y_train, k = 5),
      type.measure = "deviance",
      lambda.min.ratio = 1e-6,
      standardize = TRUE
    )
  }
)

#############################
# average curves, best lambda
#############################

lambdas = cv_runs[[1]]$lambda

# same glmnet grid in every repeat
stopifnot(all(sapply(cv_runs, function(m) identical(m$lambda, lambdas))))

cvm_mat = sapply(cv_runs, function(m) m$cvm) # rows = lambdas, cols = repeats
cvm_mean = rowMeans(cvm_mat)
cvm_sd = apply(cvm_mat, 1, sd)

i_best = which.min(cvm_mean)

# best lambda must not lie at the edge of the grid
stopifnot(i_best > 1, i_best < length(lambdas))

lambda_best = lambdas[i_best]
lambda_best

cv_curve_dat = data.frame(
  log_lambda = log(lambdas),
  cvm_mean = cvm_mean,
  cvm_sd = cvm_sd
)

plt_cv_curve = ggplot(cv_curve_dat, aes(x = log_lambda, y = cvm_mean)) +
  geom_ribbon(aes(ymin = cvm_mean - cvm_sd, ymax = cvm_mean + cvm_sd), fill = "lightgrey", alpha = 0.7) +
  geom_line(linewidth = 1.1, color = "darkblue") +
  geom_vline(xintercept = log(lambda_best), linetype = "dashed") +
  labs(
    x = "log(lambda)",
    y = paste0("Mean CV binomial deviance (", R, " x 10-fold)")
  ) +
  theme_classic()
if (interactive()) {
  plt_cv_curve
}

############################
# final model on all train
############################

fit_final = glmnet::glmnet(
  x = X_train,
  y = y_train,
  family = "binomial",
  alpha = 0,
  weights = w_train,
  lambda = lambdas,
  standardize = TRUE
)

coeffs_final = as.matrix(coef(fit_final, s = lambda_best))
coeffs_final

###################
# predict test once
###################

p_test = as.numeric(predict(fit_final, newx = X_test, s = lambda_best, type = "response"))
y_pred = ifelse(p_test > 0.5, 1, 0)

y_pred_factor = factor(ifelse(y_pred == 1, "AST", "HC"), levels = c("AST", "HC"))
y_test_factor = factor(ifelse(y_test == 1, "AST", "HC"), levels = c("AST", "HC"))

cm_test = caret::confusionMatrix(y_pred_factor, y_test_factor, positive = "AST")
cm_test

test_metrics = function(p, y) {
  pred = ifelse(p > 0.5, 1, 0)
  auc = ROCR::performance(ROCR::prediction(p, y), "auc")@y.values[[1]]
  sen = mean(pred[y == 1] == 1)
  spc = mean(pred[y == 0] == 0)
  c(auc = auc, bal_acc = (sen + spc) / 2, sen = sen, spc = spc)
}

metrics_test = test_metrics(p_test, y_test)
metrics_test

###################################
# stratified bootstrap CIs (test)
###################################

B = 2000

idx_ast = which(y_test == 1)
idx_hc = which(y_test == 0)

set.seed(202)
boot_metrics = t(sapply(
  1:B,
  function(b) {
    idx = c(
      sample(idx_ast, length(idx_ast), replace = TRUE),
      sample(idx_hc, length(idx_hc), replace = TRUE)
    )
    test_metrics(p_test[idx], y_test[idx])
  }
))

boot_ci = apply(boot_metrics, 2, quantile, probs = c(0.025, 0.975))

perf_test = data.frame(
  metric = c("AUC", "Balanced accuracy", "Sensitivity", "Specificity"),
  value = metrics_test,
  lower = boot_ci[1, ],
  upper = boot_ci[2, ]
)
rownames(perf_test) = NULL

perf_test$value_format = sapply(
  seq_along(perf_test$value),
  function(i) {
    if (perf_test$metric[i] == "AUC") {
      paste0(
        n_rnd(perf_test$value[i], 2),
        " (", n_rnd(perf_test$lower[i], 2), ", ", n_rnd(perf_test$upper[i], 2), ")"
      )
    } else {
      paste0(
        n_rnd(100 * perf_test$value[i], 1), "%",
        " (", n_rnd(100 * perf_test$lower[i], 1), ", ", n_rnd(100 * perf_test$upper[i], 1), ")"
      )
    }
  }
)
perf_test

################
# test ROC curve
################

perf_roc = ROCR::performance(ROCR::prediction(p_test, y_test), "tpr", "fpr")

roc_test_dat = data.frame(
  fpr = perf_roc@x.values[[1]],
  tpr = perf_roc@y.values[[1]]
)

roc_labels = paste0(
  "AUC: ", perf_test$value_format[1], "\n",
  "Sensitivity: ", perf_test$value_format[3], "\n",
  "Specificity: ", perf_test$value_format[4], "\n",
  "Balanced accuracy: ", perf_test$value_format[2]
)

plt_test_roc = ggplot(roc_test_dat, aes(x = fpr, y = tpr)) +
  geom_area(position = "identity", fill = "lightgrey", alpha = 0.7) +
  geom_line(data = data.frame(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), aes(x = x, y = y), linetype = "dashed") +
  geom_line(linewidth = 1.1, color = "darkblue") +
  annotate(
    "text",
    x = 0.35,
    y = 0.15,
    label = roc_labels,
    hjust = 0,
    size = 3.2
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 15, family = "sans"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )
if (interactive()) {
  plt_test_roc
}

#####################################
# predicted probabilities by group
#####################################

proba_test_dat = data.frame(
  m_id = test_dat$m_id,
  diagnosis_simple = test_dat$diagnosis_simple,
  p_test = p_test
)

plt_test_proba = ggplot(proba_test_dat, aes(x = diagnosis_simple, y = p_test)) +
  geom_boxplot(
    width = 0.5,
    fill = "white",
    linewidth = 0.4,
    alpha = 0.1,
    outlier.shape = NA
  ) +
  geom_point(aes(
    color = diagnosis_simple),
    size = 2.5,
    alpha = 0.8,
    position = position_jitter(width = 0.15, height = 0, seed = 123)
  ) +
  scale_color_manual(
    # values = c("#29ABE2", "#D9534F", "#A1A79E", "#F2D75C"),
    values = c("#CC6677", "#4477AA"),
    name = NULL
  ) +
  scale_fill_manual(
    # values = c("#29ABE2", "#D9534F", "#A1A79E", "#F2D75C"),
    values = c("#CC6677", "#4477AA"),
    name = NULL
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ylim(0, 1) +
  labs(
    x = "",
    y = "Predicted probability"
  ) +
  theme_classic()
if (interactive()) {
  plt_test_proba
}

plt_tot = plot_grid(
  plt_test_roc, plt_test_proba
)
plt_tot

################################
# permutation feature importance
################################

perm_imp = perm_importance_fn(
  model = fit_final,
  X_test = X_test,
  y_test = y_test,
  sensors = sensors,
  baseline_auc = unname(metrics_test["auc"]),
  lambda_value = lambda_best,
  B = 1000
)
perm_imp

perm_importance_plot_fn(
  perm_imp,
  range(c(perm_imp$ci_low, perm_imp$ci_high), na.rm = TRUE)
)

############
# save plots
############

# ggsave(
#   filename = "./plots/ridge_train_test_roc.png",
#   plot = plt_test_roc,
#   width = 12,
#   height = 12,
#   units = "cm",
#   dpi = 400,
#   bg = "white"
# )

########################################################
# This is AUC Shapley, NOT conventional prediction SHAP.
########################################################

# Only 13 sensors: all 2^13 = 8192 subsets can be evaluated exactly.
X_shap = X_test[, sensors, drop = FALSE]
p_shap = ncol(X_shap)
stopifnot(p_shap == 13L, all(sensors %in% rownames(coeffs_final)))

beta_shap = as.numeric(coeffs_final[sensors, 1])
names(beta_shap) = sensors

# Each column gives its contribution to the fixed model's log-odds.
# AUC is invariant to the intercept and logistic transformation.
sensor_logit = sweep(X_shap, 2, beta_shap, `*`)

masks = 0:(2^p_shap - 1L)
in_subset = outer(
  masks, 0:(p_shap - 1L),
  FUN = function(mask, bit) bitwAnd(mask, bitwShiftL(1L, bit)) != 0L
)
subset_size = rowSums(in_subset)

# AUC for every subset, evaluated on the ORIGINAL test labels.
subset_scores = in_subset %*% t(sensor_logit)
subset_auc = vapply(
  seq_len(nrow(subset_scores)),
  function(i) auc_fn(y_true = y_test, y_prob = subset_scores[i, ]),
  numeric(1)
)

# The empty subset gives tied scores (AUC = 0.5).
# The complete subset MUST reproduce the fitted ridge model's AUC.
stopifnot(abs(subset_auc[1] - 0.5) < 1e-10)
stopifnot(abs(tail(subset_auc, 1) - unname(metrics_test['auc'])) < 1e-10)

# Exact Shapley weights: 1 / (p * choose(p - 1, |S|)).
auc_shapley = vapply(seq_len(p_shap), function(j) {
  without_j = which(!in_subset[, j])
  with_j = bitwOr(masks[without_j], bitwShiftL(1L, j - 1L)) + 1L
  weights = 1 / (p_shap * choose(p_shap - 1L, subset_size[without_j]))
  sum(weights * (subset_auc[with_j] - subset_auc[without_j]))
}, numeric(1))

stopifnot(abs(sum(auc_shapley) - (unname(metrics_test['auc']) - 0.5)) < 1e-10)

auc_shapley_results = data.frame(
  sensor = sensors,
  auc_shapley = auc_shapley,
  row.names = NULL
)
auc_shapley_results = auc_shapley_results[
  order(-auc_shapley_results$auc_shapley), , drop = FALSE
]
print(auc_shapley_results)
cat('Sum of AUC Shapley values:', sum(auc_shapley), '\n')

plt_auc_shapley = ggplot(
  auc_shapley_results,
  aes(x = auc_shapley, y = reorder(sensor, auc_shapley))
) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_col(fill = 'darkblue', width = 0.7) +
  labs(x = 'Exact Shapley contribution to test AUC above 0.5', y = 'Sensor') +
  theme_classic()
print(plt_auc_shapley)

########################################################
# Standard interventional SHAP of individual predictions
# Exact for linear ridge log-odds with a training-mean baseline.
# Unlike AUC Shapley, these values explain the model's log-odds,
# NOT AUC. They do not account for conditional sensor dependence.
########################################################

reference_mean = colMeans(X_train[, sensors, drop = FALSE])
shap_logit = sweep(X_shap, 2, reference_mean, `-`)
shap_logit = sweep(shap_logit, 2, beta_shap, `*`)

baseline_logit = as.numeric(coeffs_final['(Intercept)', 1]) +
  sum(beta_shap * reference_mean)
model_logit = as.numeric(predict(
  fit_final, newx = X_test, s = lambda_best, type = 'link'
))

# Check exact reconstruction of each participant's fitted log-odds.
stopifnot(max(abs(baseline_logit + rowSums(shap_logit) - model_logit)) < 1e-8)

shap_logit_results = data.frame(
  sensor = sensors,
  mean_abs_shap_logit = colMeans(abs(shap_logit)),
  mean_shap_logit = colMeans(shap_logit),
  row.names = NULL
)
shap_logit_results = shap_logit_results[
  order(-shap_logit_results$mean_abs_shap_logit), , drop = FALSE
]
print(shap_logit_results)

plt_shap_logit = ggplot(
  shap_logit_results,
  aes(x = mean_abs_shap_logit,
      y = reorder(sensor, mean_abs_shap_logit))
) +
  geom_col(fill = 'darkblue', width = 0.7) +
  labs(x = 'Mean absolute SHAP value (log-odds)', y = 'Sensor') +
  theme_classic()
print(plt_shap_logit)
