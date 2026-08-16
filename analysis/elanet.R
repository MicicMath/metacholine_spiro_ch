source("./data_prep.R")

dat = final_ast_hc_cf_dat
dat = dat[dat$diagnosis != "CF", , drop = FALSE]

tm = Sys.time()

###############################################################
# restrict to mid and high metacholine response only for asthma
###############################################################

mid_high_metach_response_mask = dat$metacholine_response %in% c("Mid", "High")
other_mask = dat$diagnosis_simple != "AST"

dat = dat[mid_high_metach_response_mask | other_mask, , drop = FALSE]

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

X = as.matrix(dat[, sensors])
y = ifelse(dat$diagnosis_simple == "AST", 1, 0)

set.seed(101)
bundle = train_test_split(X = X, y = y, K = 5, R = 100)

ALPHAS = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
LAMBDA = "lambda.1se"
LAMBDA_IDX = ifelse(LAMBDA == "lambda.1se", "1se", "min")

########
# elanet
########

.make_inner_foldid = function(y) {
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

results_elastic_net = mclapply(
  bundle,
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
    
    inner_foldid = .make_inner_foldid(y_train)

    cv_models = lapply(
      ALPHAS,
      function(alpha) {
        inner_cv = cv.glmnet(
          x = X_train,
          y = y_train,
          foldid = inner_foldid,
          weights = w_train,
          family = "binomial",
          alpha = alpha,
          type.measure = "auc",
          standardize = TRUE
        )
      }
    )

    alpha_scores = sapply(
      cv_models,
      function(m) m$cvm[m$index[LAMBDA_IDX, 1]]
    )
    best_idx = which.max(alpha_scores)

    best_model = cv_models[[best_idx]]
    best_alpha = ALPHAS[best_idx]

    y_prob = predict(best_model, newx = X_test, s = LAMBDA, type = "response")

    y_pred_factor = factor(ifelse(y_prob > 0.5, "AST", "HC"), levels = c("AST", "HC"))
    y_test_factor = factor(ifelse(y_test == 1, "AST", "HC"), levels = c("AST", "HC"))

    cm = caret::confusionMatrix(y_pred_factor, y_test_factor, positive = "AST")

    coeffs = coef(best_model, s = LAMBDA)
    matrix_coefs = as.matrix(coeffs)
    selected_features = rownames(matrix_coefs)[matrix_coefs[, 1] != 0]

    selected_features = setdiff(selected_features, "(Intercept)")

    return(
      list(
        coeffs = matrix_coefs,
        selected_features = selected_features,
        tuned_lambda = best_model[[LAMBDA]],
        tuned_alpha = best_alpha,
        y_prob = y_prob,
        y_pred = y_pred_factor,
        y_test = y_test,
        y_test_factor = y_test_factor,
        cm = cm
      )
    )
  },
  mc.cores = parallel::detectCores() - 1
)

tm = Sys.time() - tm

#############################
# tuning porameteres selected
#############################

tuned_lambdas = sapply(results_elastic_net, function(x) x$tuned_lambda)
summary(tuned_lambdas)
tuned_alphas = sapply(results_elastic_net, function(x) x$tuned_alpha)
table(tuned_alphas)

#####################
# results performance
#####################

auc_distribution = sapply(
  results_elastic_net,
  function(x) {
    pred = ROCR::prediction(as.numeric(x$y_prob), x$y_test)
    ROCR::performance(pred, "auc")@y.values[[1]]
  }
)

auc_estimate = mean(auc_distribution)

acc_distribution = sapply(
  results_elastic_net,
  function(x) {
    x$cm$byClass[["Balanced Accuracy"]]
  }
)

acc_estimate = mean(acc_distribution)

sen_distribution = sapply(
  results_elastic_net,
  function(x) {
    x$cm$byClass[["Sensitivity"]]
  }
)

sen_estimate = mean(sen_distribution)

spc_distribution = sapply(
  results_elastic_net,
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

#############################
# feature selection frequency
#############################

selected_features = lapply(
  results_elastic_net,
  function(x) {
    x$selected_features
  }
)
selected_features = unlist(selected_features)

frequency_selected_features = as.data.frame(table(selected_features))
frequency_selected_features = frequency_selected_features[
  order(frequency_selected_features$Freq, decreasing = TRUE), ,
  drop = FALSE
]
frequency_selected_features$pct = 100 * frequency_selected_features$Freq / length(results_elastic_net)
frequency_selected_features

plt_sel_feat = ggplot(
  frequency_selected_features,
  aes(
    x = reorder(selected_features, pct),
    y = pct
  )
) +
  geom_col(
    fill = "grey55",
    width = 0.6
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct)),
    hjust = -0.1,
    size = 3.2,
    color = "grey20"
  ) +
  coord_flip() +
  labs(
    x = "",
    y = "Selection frequency"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 110),
    breaks = seq(0, 100, by = 20)
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 9),
    plot.subtitle = element_text(hjust = 9),
    text = element_text(size = 15, family = "sans"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
plt_sel_feat



