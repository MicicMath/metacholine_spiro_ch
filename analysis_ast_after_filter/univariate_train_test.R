source("./data_prep.R")

train_dat = final_ast_hc_cf_dat
test_dat = test_set_final_ast_hc_cf_dat

#######################################
# Welch comparisons within each dataset
#######################################

welch_sensor_tests = function(data) {
  results = lapply(
    sensors,
    function(x) {
      ast = data[data$diagnosis_simple == "AST", x]
      hc = data[data$diagnosis_simple == "HC", x]

      # Use available observations separately for each sensor.
      ast = ast[is.finite(ast)]
      hc = hc[is.finite(hc)]

      if (length(ast) < 2 || length(hc) < 2) {
        stop(paste("At least two observations per group are required for", x))
      }

      tt = t.test(
        ast,
        hc,
        conf.level = 0.95,
        var.equal = FALSE
      )

      # Match the unpooled-SD effect size in univariate_comparison.R.
      d_eff = effectsize::cohens_d(
        ast,
        hc,
        pooled_sd = FALSE,
        ci = 0.95
      )

      data.frame(
        sensor = x,
        n_ast = length(ast),
        n_hc = length(hc),
        mean_ast = mean(ast),
        sd_ast = sd(ast),
        mean_hc = mean(hc),
        sd_hc = sd(hc),
        mean_diff = unname(tt$estimate[[1]] - tt$estimate[[2]]),
        ci_low = tt$conf.int[[1]],
        ci_high = tt$conf.int[[2]],
        cohens_d = d_eff$Cohens_d,
        d_ci_low = d_eff$CI_low,
        d_ci_high = d_eff$CI_high,
        p_val = tt$p.value,
        row.names = NULL
      )
    }
  )

  results = do.call(rbind, results)
  # BH/FDR adjustment across sensors, separately for train and test.
  results$p_adj = p.adjust(results$p_val, method = "fdr")
  return(results)
}

t_test_train = welch_sensor_tests(train_dat)
t_test_test = welch_sensor_tests(test_dat)

##################
# formatted tables
##################

format_welch_results = function(results) {
  data.frame(
    Sensor = results$sensor,
    "AST n" = results$n_ast,
    "HC n" = results$n_hc,
    "AST mean (SD)" = paste0(
      n_rnd(results$mean_ast, 2), " (", n_rnd(results$sd_ast, 2), ")"
    ),
    "HC mean (SD)" = paste0(
      n_rnd(results$mean_hc, 2), " (", n_rnd(results$sd_hc, 2), ")"
    ),
    "AST-HC mean difference (95% CI)" = paste0(
      n_rnd(results$mean_diff, 2), " (",
      n_rnd(results$ci_low, 2), ", ", n_rnd(results$ci_high, 2), ")"
    ),
    "Cohen's d (95% CI)" = paste0(
      n_rnd(results$cohens_d, 2), " (",
      n_rnd(results$d_ci_low, 2), ", ", n_rnd(results$d_ci_high, 2), ")"
    ),
    "Raw p-value" = p_rnd(results$p_val),
    "FDR-adjusted p-value" = p_rnd(results$p_adj),
    check.names = FALSE,
    row.names = NULL
  )
}

p_dat_train_format = format_welch_results(t_test_train)
p_dat_test_format = format_welch_results(t_test_test)

cat("\nTraining: AST vs HC (Welch)\n")
print(p_dat_train_format, row.names = FALSE)

cat("\nTest: AST vs HC (Welch)\n")
print(p_dat_test_format, row.names = FALSE)

#####################################
# boxplots: paired train/test panels
#####################################

plot_cols = c("m_id", "diagnosis_simple", sensors)

plot_dat = rbind(
  data.frame(dataset = "Training", train_dat[, plot_cols]),
  data.frame(dataset = "Test", test_dat[, plot_cols])
)
plot_dat$dataset = factor(plot_dat$dataset, levels = c("Training", "Test"))
plot_dat$diagnosis_simple = factor(plot_dat$diagnosis_simple, levels = c("AST", "HC"))

t_test_all = rbind(
  data.frame(dataset = "Training", t_test_train),
  data.frame(dataset = "Test", t_test_test)
)

diagnosis_colors = c(AST = "#CC6677", HC = "#4477AA")

bps = lapply(
  sensors,
  function(x) {
    sensor_dat = plot_dat[is.finite(plot_dat[[x]]), , drop = FALSE]
    sensor_tests = t_test_all[t_test_all$sensor == x, , drop = FALSE]

    # Shared scale and bracket height across train/test for this sensor.
    y_range = range(sensor_dat[[x]])
    y_position = y_range[2] + 0.08 * diff(y_range)

    pwc = data.frame(
      dataset = factor(sensor_tests$dataset, levels = c("Training", "Test")),
      group1 = "AST",
      group2 = "HC",
      y.position = y_position,
      p.adj = sensor_tests$p_adj,
      p.adj.str = paste0(
        # "FDR p ",
        ifelse(sensor_tests$p_adj < 0.001, "", "= "),
        p_rnd(sensor_tests$p_adj)
      )
    )

    bp = ggboxplot(sensor_dat,
      x = "diagnosis_simple", y = x,
      width = 0.5,
      fill = "white",
      linewidth = 0.4,
      alpha = 0.1, legend = "none", outlier.shape = NA
    ) +
      geom_point(
        aes(x = diagnosis_simple, y = !!sym(x), fill = diagnosis_simple, color = diagnosis_simple),
        size = 2.4,
        alpha = 0.5,
        position = position_jitter(width = 0.15, height = 0, seed = 123)
      ) +
      scale_color_manual(values = diagnosis_colors, name = NULL) +
      scale_fill_manual(values = diagnosis_colors, name = NULL) +
      facet_wrap(~dataset, nrow = 1, scales = "fixed") +
      stat_pvalue_manual(
        pwc,
        label = "p.adj.str",
        tip.length = 0.01,
        label.size = 3,
        hide.ns = FALSE
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(
        subtitle = paste0("Sensor: ", x),
        x = "",
        y = "Intensity [a.u.]"
      ) +
      theme_classic() +
      ggtheme_no_legend +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(size = 10)
      )

    bp
  }
)
names(bps) = sensors

plt_bps = plot_grid(
  plotlist = bps,
  ncol = 5,
  align = "hv"
)
print(plt_bps)

###############
# optional saves
###############

# dir.create("./tables", showWarnings = FALSE)
# write.csv(p_dat_train_format, "./tables/univariate_train_welch.csv", row.names = FALSE)
# write.csv(p_dat_test_format, "./tables/univariate_test_welch.csv", row.names = FALSE)

# ggsave(
#   filename = "./plots/univariate_train_test_boxplots.png",
#   plot = plt_bps,
#   width = 36,
#   height = 42,
#   units = "cm",
#   dpi = 300,
#   bg = "white"
# )
