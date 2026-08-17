source("./data_prep.R")

dat = final_ast_hc_cf_dat
dat = dat[dat$diagnosis != "CF", , drop = FALSE]

###############################################################
# restrict to mid and high metacholine response only for asthma
###############################################################

mid_high_metach_response_mask = dat$metacholine_response %in% c("Low", "Mid", "High")
other_mask = dat$diagnosis_simple != "AST"

dat = dat[mid_high_metach_response_mask | other_mask, , drop = FALSE]

table(dat$diagnosis_simple)

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

X_dat = dat[, sensors, drop = FALSE]
y_dat = ifelse(dat$diagnosis_simple == "AST", 1, 0)

set.seed(101)
bundle = train_test_split(X = X, y = y, K = 5, R = 100)

#######
# svm
#######

X = as.matrix(dat[, sensors])
y = ifelse(dat$diagnosis_simple == "AST", 1, 0)

set.seed(101)
bundle = train_test_split(X = X, y = y, K = 5, R = 100)

results_svm = mclapply(
  bundle,
  function(split) {
    X_train = split$X_train
    y_train = split$y_train
    X_test = split$X_test
    y_test = split$y_test

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

    cost_grid = 2^seq(-5, 3, by = 1)

    tuning_auc = sapply(
      cost_grid,
      function(j) {
        auc_folds = sapply(
          sort(unique(inner_foldid)),
          function(k) {
            train_idx = inner_foldid != k
            test_idx = inner_foldid == k

            y_inner_train = y_train[train_idx]

            y_inner_train_factor = factor(
              ifelse(y_inner_train == 1, "AST", "HC"),
              levels = c("AST", "HC")
            )

            class_weights = c(
              AST = 0.5 / sum(y_inner_train == 1),
              HC = 0.5 / sum(y_inner_train == 0)
            )

            class_weights = class_weights / mean(class_weights)

            model = e1071::svm(
              x = X_train[train_idx, , drop = FALSE],
              y = y_inner_train_factor,
              type = "C-classification",
              kernel = "linear",
              cost = j,
              class.weights = class_weights,
              # scale = TRUE,
              probability = TRUE
            )

            pred = predict(
              model,
              X_train[test_idx, , drop = FALSE],
              probability = TRUE
            )

            probas = attr(pred, "probabilities")[, "AST"]

            pred_auc = ROCR::prediction(
              probas,
              y_train[test_idx]
            )

            auc = ROCR::performance(
              pred_auc,
              "auc"
            )@y.values[[1]]

            return(auc)
          }
        )

        return(mean(auc_folds))
      }
    )

    best_idx = which.max(tuning_auc)

    best_cost = cost_grid[best_idx]

    y_train_factor = factor(
      ifelse(y_train == 1, "AST", "HC"),
      levels = c("AST", "HC")
    )

    class_weights = c(
      AST = 0.5 / sum(y_train == 1),
      HC = 0.5 / sum(y_train == 0)
    )

    class_weights = class_weights / mean(class_weights)

    tuned_model = e1071::svm(
      x = X_train,
      y = y_train_factor,
      type = "C-classification",
      kernel = "linear",
      cost = best_cost,
      class.weights = class_weights,
      # scale = TRUE,
      probability = TRUE
    )

    pred_class = predict(
      tuned_model,
      X_test,
      probability = TRUE
    )

    probas = attr(pred_class, "probabilities")[, "AST"]

    predicted = ifelse(pred_class == "AST", 1, 0)

    y_pred_factor = factor(
      ifelse(predicted == 1, "AST", "HC"),
      levels = c("AST", "HC")
    )

    y_test_factor = factor(
      ifelse(y_test == 1, "AST", "HC"),
      levels = c("AST", "HC")
    )

    cm = caret::confusionMatrix(
      y_pred_factor,
      y_test_factor,
      positive = "AST"
    )

    return(
      list(
        cm = cm,
        y_test = y_test,
        y_proba = probas,
        y_pred = predicted,
        best_cost = best_cost
      )
    )
  }, mc.cores = parallel::detectCores() - 2
)

#####################
# results performance
#####################

auc_distribution = sapply(
  results_svm,
  function(x) {
    pred = ROCR::prediction(
      as.numeric(x$y_proba),
      x$y_test
    )

    ROCR::performance(pred, "auc")@y.values[[1]]
  }
)

auc_estimate = mean(auc_distribution)

acc_distribution = sapply(
  results_svm,
  function(x) {
    x$cm$byClass[["Balanced Accuracy"]]
  }
)

acc_estimate = mean(acc_distribution)

sen_distribution = sapply(
  results_svm,
  function(x) {
    x$cm$byClass[["Sensitivity"]]
  }
)

sen_estimate = mean(sen_distribution)

spc_distribution = sapply(
  results_svm,
  function(x) {
    x$cm$byClass[["Specificity"]]
  }
)

spc_estimate = mean(spc_distribution)

perf_dat = data.frame(
  metric = c(
    "AUC",
    "Balanced accuracy",
    "Sensitivity",
    "Specificity"
  ),
  value = c(
    auc_estimate,
    100 * acc_estimate,
    100 * sen_estimate,
    100 * spc_estimate
  ),
  sds = c(
    sd(auc_distribution),
    sd(acc_distribution * 100),
    sd(sen_distribution * 100),
    sd(spc_distribution * 100)
  )
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
# tuning parameters
####################

best_cost = sapply(
  results_svm,
  function(x) x$best_cost
)

table(best_cost)

#######
# rf
#######

results_rf = mclapply(
  bundle,
  function(split) {
    X_train = split$X_train
    y_train = split$y_train
    X_test = split$X_test
    y_test = split$y_test

    # .make_foldid = function(y) {
    #   folds = caret::createFolds(
    #     factor(y, levels = c(1, 0)),
    #     k = 5,
    #     returnTrain = FALSE
    #   )

    #   foldid = rep(0, length(y))

    #   for (i in seq_along(folds)) {
    #     foldid[folds[[i]]] = i
    #   }

    #   return(foldid)
    # }

    # inner_foldid = .make_foldid(y_train)

    # num_trees_grid = c(500, 1000, 2000)

    mtry = floor(sqrt(ncol(X_train)))

    # tuning_auc = sapply(
    #   num_trees_grid,
    #   function(j) {
    #     auc_folds = sapply(
    #       sort(unique(inner_foldid)),
    #       function(k) {
    #         train_idx = inner_foldid != k
    #         test_idx = inner_foldid == k

    #         y_inner_train = y_train[train_idx]

    #         y_inner_train_factor = factor(
    #           ifelse(y_inner_train == 1, "AST", "HC"),
    #           levels = c("AST", "HC")
    #         )

    #         class_weights = c(
    #           AST = 0.5 / sum(y_inner_train == 1),
    #           HC = 0.5 / sum(y_inner_train == 0)
    #         )

    #         class_weights = class_weights / mean(class_weights)

    #         model = ranger::ranger(
    #           x = X_train[train_idx, , drop = FALSE],
    #           y = y_inner_train_factor,
    #           num.trees = j,
    #           mtry = mtry,
    #           min.node.size = 1,
    #           class.weights = class_weights,
    #           probability = TRUE,
    #           num.threads = 1
    #         )

    #         probas = predict(
    #           model,
    #           X_train[test_idx, , drop = FALSE]
    #         )$predictions[, "AST"]

    #         pred_auc = ROCR::prediction(
    #           probas,
    #           y_train[test_idx]
    #         )

    #         auc = ROCR::performance(
    #           pred_auc,
    #           "auc"
    #         )@y.values[[1]]

    #         return(auc)
    #       }
    #     )

    #     return(mean(auc_folds))
    #   }
    # )

    # best_idx = which.max(tuning_auc)

    # best_num_trees = num_trees_grid[best_idx]

    y_train_factor = factor(
      ifelse(y_train == 1, "AST", "HC"),
      levels = c("AST", "HC")
    )

    class_weights = c(
      AST = 0.5 / sum(y_train == 1),
      HC = 0.5 / sum(y_train == 0)
    )

    class_weights = class_weights / mean(class_weights)

    tuned_model = ranger::ranger(
      x = X_train,
      y = y_train_factor,
      num.trees = 500,
      mtry = mtry,
      min.node.size = 1,
      class.weights = class_weights,
      probability = TRUE,
      num.threads = 1
    )

    probas = predict(
      tuned_model,
      X_test
    )$predictions[, "AST"]

    predicted = ifelse(probas > 0.5, 1, 0)

    y_pred_factor = factor(
      ifelse(predicted == 1, "AST", "HC"),
      levels = c("AST", "HC")
    )

    y_test_factor = factor(
      ifelse(y_test == 1, "AST", "HC"),
      levels = c("AST", "HC")
    )

    cm = caret::confusionMatrix(
      y_pred_factor,
      y_test_factor,
      positive = "AST"
    )

    return(
      list(
        cm = cm,
        y_test = y_test,
        y_proba = probas,
        y_pred = predicted
        # best_num_trees = best_num_trees
      )
    )
  } , mc.cores = parallel::detectCores() - 2
)

#####################
# results performance
#####################

auc_distribution = sapply(
  results_rf,
  function(x) {
    pred = ROCR::prediction(
      as.numeric(x$y_proba),
      x$y_test
    )

    ROCR::performance(pred, "auc")@y.values[[1]]
  }
)

auc_estimate = mean(auc_distribution)

acc_distribution = sapply(
  results_rf,
  function(x) {
    x$cm$byClass[["Balanced Accuracy"]]
  }
)

acc_estimate = mean(acc_distribution)

sen_distribution = sapply(
  results_rf,
  function(x) {
    x$cm$byClass[["Sensitivity"]]
  }
)

sen_estimate = mean(sen_distribution)

spc_distribution = sapply(
  results_rf,
  function(x) {
    x$cm$byClass[["Specificity"]]
  }
)

spc_estimate = mean(spc_distribution)

perf_dat = data.frame(
  metric = c(
    "AUC",
    "Balanced accuracy",
    "Sensitivity",
    "Specificity"
  ),
  value = c(
    auc_estimate,
    100 * acc_estimate,
    100 * sen_estimate,
    100 * spc_estimate
  ),
  sds = c(
    sd(auc_distribution),
    sd(acc_distribution * 100),
    sd(sen_distribution * 100),
    sd(spc_distribution * 100)
  )
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
# tuning parameters
####################

# best_num_trees = sapply(
#   results_rf,
#   function(x) x$best_num_trees
# )
# 
# table(best_num_trees)
