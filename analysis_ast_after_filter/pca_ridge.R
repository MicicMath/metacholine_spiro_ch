source("./data_prep.R")

dat = final_ast_hc_cf_dat
table(dat$diagnosis_simple)

#####################
# restrict to batch 2
#####################

# dat = dat[dat$batch == "batch_20260808", , drop = FALSE]

# subset_sensors = c(
#   # "S1",
#   # "S3",
#   # "S4",
#   # "S5",
#   # "S6",
#   # "S7",
#   "S1BH",
#   "S2BH",
#   "S3BH",
#   "S4BH",
#   "S5BH",
#   "S6BH",
#   "S7BH"
# )
# 
# sensors = setdiff(sensors, subset_sensors)
# 
# cols = setdiff(colnames(dat), subset_sensors)
# 
# dat = dat[, cols]

#############
# pca ridge
#############

train_test_split = function(X, y, K = 5, R = 100) {
  flds = lapply(
    1:R,
    function(i) {
      createFolds(factor(y, levels = c(1, 0)), k = K, returnTrain = FALSE)
    }
  )
  flds = unlist(flds, recursive = FALSE)
  names(flds) = paste0("f", 1:length(flds))

  splits_bundle = lapply(
    flds,
    function(x) {
      X_train = X[-x, , drop = FALSE]
      X_test = X[x, , drop = FALSE]

      y_train = y[-x]
      y_test = y[x]

      return(list(
        X_train = X_train,
        X_test = X_test,
        y_train = y_train,
        y_test = y_test
      ))
    }
  )
  return(splits_bundle)
}

X = as.matrix(dat[, sensors])
y = ifelse(dat$diagnosis_simple == "AST", 1, 0)

set.seed(101)
bundle = train_test_split(X = X, y = y, K = 5, R = 100)

lambda_grid = 10^seq(-4, 4, length.out = 100)

results_pca_ridge = mclapply(
  seq_along(bundle),
  function(i) {
    set.seed(1000 + i)

    split = bundle[[i]]

    X_train = split$X_train
    y_train = split$y_train
    X_test = split$X_test
    y_test = split$y_test

    r = 0.5 # 0.5 is standard balance
    w_train = ifelse(y_train == 1,
      r * (1 / sum(y_train == 1)),
      (1 - r) * (1 / sum(y_train == 0))
    )

    w_train = w_train * (length(y_train) / sum(w_train))

    y_train_factor = factor(ifelse(y_train == 1, "AST", "HC"), levels = c("AST", "HC"))

    inner_index = caret::createFolds(
      y_train_factor,
      k = 5,
      returnTrain = TRUE
    )

    one_se_lambda = function(x, metric, maximize) {
      best_idx = if (maximize) {
        which.max(x[[metric]])
      } else {
        which.min(x[[metric]])
      }

      best_metric = x[[metric]][best_idx]
      metric_se = x[[paste0(metric, "SD")]][best_idx] / sqrt(length(inner_index))

      if (is.na(metric_se)) {
        return(best_idx)
      }

      if (maximize) {
        metric_cutoff = best_metric - metric_se
        eligible = which(x[[metric]] >= metric_cutoff)
      } else {
        metric_cutoff = best_metric + metric_se
        eligible = which(x[[metric]] <= metric_cutoff)
      }

      eligible[which.max(x$lambda[eligible])]
    }

    ctrl = caret::trainControl(
      method = "cv",
      number = 5,
      index = inner_index,
      classProbs = TRUE,
      summaryFunction = caret::twoClassSummary,
      selectionFunction = one_se_lambda,
      savePredictions = "final",
      allowParallel = FALSE,
      preProcOptions = list(thresh = 0.8)
    )

    tuned_model = caret::train(
      x = X_train,
      y = y_train_factor,
      method = "glmnet",
      family = "binomial",
      metric = "ROC",
      maximize = TRUE,
      weights = w_train,
      trControl = ctrl,
      preProcess = c("center", "scale", "pca"),
      tuneGrid = expand.grid(
        alpha = 0,
        lambda = lambda_grid
      )
    )

    probas = predict(tuned_model, newdata = X_test, type = "prob")[, "AST"]
    predicted = ifelse(probas > 0.5, 1, 0)

    n_pc = if (!is.null(tuned_model$preProcess$numComp)) {
      tuned_model$preProcess$numComp
    } else {
      ncol(predict(tuned_model$preProcess, X_train[1, , drop = FALSE]))
    }

    y_pred_factor = factor(ifelse(predicted == 1, "AST", "HC"), levels = c("AST", "HC"))
    y_test_factor = factor(ifelse(y_test == 1, "AST", "HC"), levels = c("AST", "HC"))

    cm = caret::confusionMatrix(y_pred_factor, y_test_factor, positive = "AST")

    return(
      list(
        cm = cm,
        tuned_lambda = tuned_model$bestTune$lambda,
        n_pc = n_pc,
        y_test = y_test,
        y_proba = probas,
        y_pred = predicted
      )
    )
  }, mc.cores = parallel::detectCores() - 2
)

###################
# tuned parameters
###################

tuned_lambdas = sapply(
  results_pca_ridge,
  function(x) {
    x$tuned_lambda
  }
)

summary(tuned_lambdas)

n_pc = sapply(
  results_pca_ridge,
  function(x) {
    x$n_pc
  }
)

table(n_pc)
summary(n_pc)

#####################
# results performance
#####################

auc_distribution = sapply(
  results_pca_ridge,
  function(x) {
    pred = ROCR::prediction(as.numeric(x$y_proba), x$y_test)
    ROCR::performance(pred, "auc")@y.values[[1]]
  }
)

auc_estimate = mean(auc_distribution)

acc_distribution = sapply(
  results_pca_ridge,
  function(x) {
    x$cm$byClass[["Balanced Accuracy"]]
  }
)

acc_estimate = mean(acc_distribution)

sen_distribution = sapply(
  results_pca_ridge,
  function(x) {
    x$cm$byClass[["Sensitivity"]]
  }
)

sen_estimate = mean(sen_distribution)

spc_distribution = sapply(
  results_pca_ridge,
  function(x) {
    x$cm$byClass[["Specificity"]]
  }
)

spc_estimate = mean(spc_distribution)

perf_dat = data.frame(
  metric = c("AUC", "Balanced accuracy", "Sensitivity", "Specificity"),
  value = c(auc_estimate, 100 * acc_estimate, 100 * sen_estimate, 100 * spc_estimate),
  sds = c(sd(auc_distribution), sd(acc_distribution * 100), sd(sen_distribution * 100), sd(spc_distribution * 100))
)

perf_dat$value_format = sapply(
  seq_along(perf_dat$value),
  function(i) {
    x = perf_dat$value[i]
    y = perf_dat$sds[i]
    name = perf_dat$metric[i]
    if (name == "AUC") {
      res = sprintf("%.2f", x)
      std = sprintf("%.3f", y)
      return(paste0(res, " (SD ", std, ")"))
    } else {
      res = paste0(sprintf("%.1f", x), "%")
      std = sprintf("%.1f", y)
      return(paste0(res, " (SD ", std, ")"))
    }
  }
)
perf_dat

####################
# average ROC curve
####################

fpr_grid = seq(0, 1, by = 0.05)

roc_tpr = sapply(
  seq_along(results_pca_ridge),
  function(i) {
    pred = ROCR::prediction(
      as.numeric(results_pca_ridge[[i]]$y_proba),
      results_pca_ridge[[i]]$y_test
    )

    perf = ROCR::performance(pred, "tpr", "fpr")

    roc_dat = data.frame(
      fpr = perf@x.values[[1]],
      tpr = perf@y.values[[1]]
    )

    approx(
      x = roc_dat$fpr,
      y = roc_dat$tpr,
      xout = fpr_grid,
      ties = max
    )$y
  }
)

mean_roc_dat = data.frame(
  fpr = fpr_grid,
  tpr = rowMeans(roc_tpr, na.rm = TRUE)
)

mean_roc_dat = rbind(
  data.frame(fpr = 0, tpr = 0),
  mean_roc_dat,
  data.frame(fpr = 1, tpr = 1)
)

roc_labels = paste0(
  "AUC: ", n_rnd(auc_estimate, n = 2), "\n",
  "Sensitivity: ", n_rnd(sen_estimate * 100, n = 1), "%\n",
  "Specificity: ", n_rnd(spc_estimate * 100, n = 1), "%\n",
  "Balanced accuracy: ", n_rnd(acc_estimate * 100, n = 1), "%"
)

plt_mean_roc = ggplot(mean_roc_dat, aes(x = fpr, y = tpr)) +
  geom_area(position = "identity", fill = "lightgrey", alpha = 0.7) +
  geom_line(data = data.frame(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)), aes(x = x, y = y), linetype = "dashed") +
  geom_line(linewidth = 1.1, color = "darkblue") +
  annotate(
    "text",
    x = 0.45,
    y = 0.15,
    label = roc_labels,
    hjust = 0,
    # vjust = -0.8,
    size = 3.6
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
    plot.title = element_text(hjust = 0.5, size = 9),
    plot.subtitle = element_text(hjust = 9),
    text = element_text(size = 15, family = "sans"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )

if (interactive()) {
  plt_mean_roc
}

###########
# save plot
###########

# ggsave(
#   filename = "./plots/pca_ridge_results.png",
#   plot = plt_mean_roc,
#   width = 12,
#   height = 12,
#   units = "cm",
#   dpi = 400,
#   bg = "white"
# )
