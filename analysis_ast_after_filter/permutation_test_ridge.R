source("./data_prep.R")

tm = Sys.time()

dat = final_ast_hc_cf_dat
table(dat$diagnosis_simple)

if (!exists("observed_R")) {
  observed_R = 100
}

if (!exists("n_permutations")) {
  n_permutations = 1000
}

if (!exists("permutation_R")) {
  permutation_R = 1
}

#######
# setup
#######

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

fit_ridge_split = function(split, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

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
    alpha = 0,
    standardize = TRUE
  )

  probas = predict(tuned_model, newx = X_test, s = "lambda.1se", type = "response")
  predicted = ifelse(probas > 0.5, 1, 0)

  y_pred_factor = factor(ifelse(predicted == 1, "AST", "HC"), levels = c("AST", "HC"))
  y_test_factor = factor(ifelse(y_test == 1, "AST", "HC"), levels = c("AST", "HC"))

  cm = caret::confusionMatrix(y_pred_factor, y_test_factor, positive = "AST")

  return(
    list(
      cm = cm,
      tuned_lambda = tuned_model$lambda.1se,
      y_test = y_test,
      y_proba = probas,
      y_pred = predicted
    )
  )
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

auc_from_results = function(per_fold_results) {
  auc_distribution = sapply(
    per_fold_results,
    function(x) {
      pred = ROCR::prediction(as.numeric(x$y_proba), x$y_test)
      ROCR::performance(pred, "auc")@y.values[[1]]
    }
  )

  mean(auc_distribution)
}

X = as.matrix(dat[, sensors])
y = ifelse(dat$diagnosis_simple == "AST", 1, 0)

################
# observed ridge
################

set.seed(101)
bundle_observed = train_test_split(X = X, y = y, K = 5, R = observed_R)

results_ridge_observed = mclapply(
  seq_along(bundle_observed),
  function(i) {
    fit_ridge_split(bundle_observed[[i]], seed = 1000 + i)
  }, mc.cores = parallel::detectCores() - 2
)

perf_ridge_observed = perf_report(results_ridge_observed, "Ridge observed")
perf_ridge_observed

auc_observed = auc_from_results(results_ridge_observed)
auc_observed

####################
# permutation testing
####################

permutation_results = mclapply(
  1:n_permutations,
  function(i) {
    set.seed(2000 + i)

    y_perm = sample(y)

    set.seed(3000 + i)
    bundle_perm = train_test_split(X = X, y = y_perm, K = 5, R = permutation_R)

    results_perm = lapply(
      seq_along(bundle_perm),
      function(j) {
        fit_ridge_split(bundle_perm[[j]], seed = 400000 + i * 100 + j)
      }
    )

    auc_perm = auc_from_results(results_perm)

    return(list(
      auc = auc_perm,
      results = results_perm
    ))
  }, mc.cores = parallel::detectCores() - 2
)

auc_permutation_distribution = sapply(
  permutation_results,
  function(x) {
    x$auc
  }
)

summary(auc_permutation_distribution)

p_auc = (sum(auc_permutation_distribution >= auc_observed) + 1) /
  (length(auc_permutation_distribution) + 1)

p_auc

tm = Sys.time() - tm
