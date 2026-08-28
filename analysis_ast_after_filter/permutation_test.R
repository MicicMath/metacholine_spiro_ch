source("./data_prep.R")

tm = Sys.time()

dat = final_ast_hc_cf_dat
table(dat$diagnosis_simple)


#######
# setup
#######

train_test_split_perm = function(X, y, R = 1000, p = 0.8) {

  splits_bundle = lapply(
    1:R,
    function(i) {

      y_perm = sample(y)

      train_idx = caret::createDataPartition(
        factor(y_perm, levels = c(1, 0)),
        p = p,
        list = FALSE
      )

      list(
        X_train = X[train_idx, , drop = FALSE],
        X_test = X[-train_idx, , drop = FALSE],
        y_train = y_perm[train_idx],
        y_test = y_perm[-train_idx]
      )
    }
  )

  return(splits_bundle)
}

perf_report = function(per_fold_results, model_name) {
  auc_distribution = sapply(
    per_fold_results,
    function(x) {
      pred = ROCR::prediction(as.numeric(x$y_proba), x$y_test)
      ROCR::performance(pred, "auc")@y.values[[1]]
    }
  )

  auc_estimate = mean(auc_distribution)

  acc_distribution = sapply(
    per_fold_results,
    function(x) {
      x$cm$byClass[["Balanced Accuracy"]]
    }
  )

  acc_estimate = mean(acc_distribution)

  sen_distribution = sapply(
    per_fold_results,
    function(x) {
      x$cm$byClass[["Sensitivity"]]
    }
  )

  sen_estimate = mean(sen_distribution)

  spc_distribution = sapply(
    per_fold_results,
    function(x) {
      x$cm$byClass[["Specificity"]]
    }
  )

  spc_estimate = mean(spc_distribution)

  res = data.frame(
    model = model_name,
    auc = paste0(sprintf("%.2f", auc_estimate), " (SD ", sprintf("%.2f", sd(auc_distribution)), ")"),
    acc = paste0(sprintf("%.1f", acc_estimate * 100), "%", " (SD ", sprintf("%.1f", sd(acc_distribution * 100)), ")"),
    sen = paste0(sprintf("%.1f", sen_estimate * 100), "%", " (SD ", sprintf("%.1f", sd(sen_distribution * 100)), ")"),
    spc = paste0(sprintf("%.1f", spc_estimate * 100), "%", " (SD ", sprintf("%.1f", sd(spc_distribution * 100)), ")")
  )
  colnames(res) = c(
    "Model",
    "AUC",
    "Balanced accuracy",
    "Sensitivity",
    "Specificity"
  )
  res

}

X = as.matrix(dat[, sensors])
y = ifelse(dat$diagnosis_simple == "AST", 1, 0)

set.seed(101)

bundle_perm = train_test_split_perm(
  X = X,
  y = y,
  R = 1000,
  p = 0.8
)

#######
# lasso
#######

results_lasso = mclapply(
  bundle_perm,
  function(split) {
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

    .make_foldid = function(y) {
      folds = caret::createFolds(
        factor(y, levels = c(1, 0)),
        k = 5,
        returnTrain = FALSE
      )

      foldid = rep(0, length(y))

      for (i in seq_along(folds)) {
        foldid[folds[[i]]] = i
      }

      return(foldid)
    }

    inner_foldid = .make_foldid(y_train)

    tuned_model = glmnet::cv.glmnet(
      x = X_train,
      y = y_train,
      family = "binomial",
      type.measure = "auc",
      weights = w_train,
      foldid = inner_foldid,
      alpha = 1,
      standardize = TRUE
    )

    probas = predict(tuned_model, newx = X_test, s = "lambda.1se", type = "response")
    predicted = ifelse(probas > 0.5, 1, 0)
    coeffs = coef(tuned_model, s = "lambda.1se")

    matrix_coefs = as.matrix(coeffs)
    selected_features = rownames(matrix_coefs)[matrix_coefs[, 1] != 0]
    selected_features = setdiff(selected_features, "(Intercept)")

    y_pred_factor = factor(ifelse(predicted == 1, "AST", "HC"), levels = c("AST", "HC"))
    y_test_factor = factor(ifelse(y_test == 1, "AST", "HC"), levels = c("AST", "HC"))

    cm = caret::confusionMatrix(y_pred_factor, y_test_factor, positive = "AST")

    return(
      list(
        cm = cm,
        selected_features = selected_features,
        y_test = y_test,
        y_proba = probas,
        y_pred = predicted
      )
    )
  }, mc.cores = parallel::detectCores() - 2
)

#####################
# results performance
#####################

auc_distribution = sapply(
  results_lasso,
  function(x) {
    pred = ROCR::prediction(as.numeric(x$y_proba), x$y_test)
    ROCR::performance(pred, "auc")@y.values[[1]]
  }
)

auc_estimate = mean(auc_distribution)

acc_distribution = sapply(
  results_lasso,
  function(x) {
    x$cm$byClass[["Balanced Accuracy"]]
  }
)

acc_estimate = mean(acc_distribution)

sen_distribution = sapply(
  results_lasso,
  function(x) {
    x$cm$byClass[["Sensitivity"]]
  }
)

sen_estimate = mean(sen_distribution)

spc_distribution = sapply(
  results_lasso,
  function(x) {
    x$cm$byClass[["Specificity"]]
  }
)

spc_estimate = mean(spc_distribution)

perf_dat = data.frame(
  metric = c("AUC", "Balanced accuracy", "Sensitivity", "Specificity"),
  value = c(auc_estimate, 100 * acc_estimate, 100 * sen_estimate, 100 * spc_estimate),
  sds = c(sd(auc_distribution), sd(acc_distribution * 100), sd(sen_distribution * 100), sd(spc_distribution * 100)),
  algorithm = "LASSO"
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

perf_dat_lasso = perf_dat

auc_observed = 0.95

p_auc = (sum(auc_distribution >= auc_observed) + 1) /
  (length(auc_distribution) + 1)

p_auc

tm = Sys.time() - tm
