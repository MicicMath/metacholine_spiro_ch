source("./data_prep.R")

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
      foldid = make_foldid(y_train, k = 10),
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
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 123)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ylim(0, 1) +
  labs(
    x = "",
    y = "Predicted probability AST (test)"
  ) +
  theme_classic()
if (interactive()) {
  plt_test_proba
}

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
