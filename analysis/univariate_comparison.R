source("./data_prep.R")

dat = final_ast_hc_cf_dat[, cols]
dat = dat[dat$diagnosis != "PUL_CF", , drop = FALSE]

dat$diagnosis = ifelse(
  dat$diagnosis == "PUL_ASTHMA_ALLERGIC",
  "AST_AL",
  ifelse(dat$diagnosis == "PUL_ASTHMA_NONALLERGIC",
    "AST_NOAL", "HC"
  )
)

#############################
# exclude asthma non-allergic
#############################

dat = dat[dat$diagnosis != "AST_NOAL", , drop = FALSE]

#######################
# remap diagnosis label
#######################

dat$diagnosis_simple = ifelse(dat$diagnosis == "HC", "HC", "AST")

#######################
# add metacholine label
#######################

dat$metacho = ifelse(
  (dat$diagnosis_simple == "AST") & (dat$metacholine_test_result == 0),
  "R0", ifelse(
    (dat$diagnosis_simple == "AST") & (dat$metacholine_test_result == 1),
    "R1", ifelse(
      (dat$diagnosis_simple == "AST") & (dat$metacholine_test_result == 2),
      "R2", ifelse(
        (dat$diagnosis_simple == "AST") & (dat$metacholine_test_result == 3),
        "R3", "HC"
      )
    )
  )
)

###################################
# remove asthma with no metacholine
###################################

mask_no_metacho = (dat$diagnosis_simple == "AST") & is.na(dat$metacholine_test_result)

dat = dat[!mask_no_metacho, , drop = FALSE]

#################################
# restrict to confirmed diagnosis
#################################

conf_diag = (dat$diagnosis_status == 1) | is.na(dat$diagnosis_status)
dat = dat[conf_diag, , drop = FALSE]

#######################################
# restrict to metacholine level 2 and 3
#######################################

metacho_23 = (dat$metacho %in% c("R2", "R3")) | dat$diagnosis == "HC"
dat = dat[metacho_23, , drop = FALSE]

####
# st
####

w_test = mclapply(
  sensors,
  function(x) {
    ast = dat[dat$diagnosis_simple == "AST", x]
    hc = dat[dat$diagnosis_simple == "HC", x]
    w = wilcox.test(ast, hc, conf.int = TRUE)
    d = data.frame(
      sensor = x,
      diff_cf_vs_hc_ci = paste0(
        n_rnd(w$estimate, 2),
        " (",
        n_rnd(w$conf.int[[1]], 2),
        ", ",
        n_rnd(w$conf.int[[2]], 2),
        ")"
      ),
      p_val = w$p.value
    )
    rownames(d) = NULL
    d
  }, mc.cores = parallel::detectCores() - 2
)
w_test = do.call(rbind, w_test)

p_adj = p.adjust(w_test$p_val, method = "fdr")
names(p_adj) = sensors
w_test$p_adj = p_adj

p_dat_format = w_test
p_dat_format$p_val = p_rnd(p_dat_format$p_val)
p_dat_format$p_adj = p_rnd(p_dat_format$p_adj)

# p_dat_format[
#   ,
#   c(
#     "sensor",
#     "diff_cf_vs_hc_ci_base",
#     "p_adj_base",
#     "diff_cf_vs_hc_ci_end",
#     "p_adj_end"
#   )
# ]

colnames(p_dat_format) = c(
  "Sensor",
  "Est. difference (95% CI)",
  "Raw p-value",
  "Adjusted p-value"
)

####
# bp
####

bps = lapply(
  sensors,
  function(x) {
    pwc = data.frame(
      group1 = "AST",
      group2 = "HC",
      y.position = max(dat[[x]], na.rm = TRUE) * 1.1,
      p.adj = as.numeric(p_adj[x]),
      p.adj.str = p_rnd(p_adj[x])
    )
    bp = ggboxplot(dat,
      x = "diagnosis_simple", y = x,
      width = 0.5,
      fill = "white",
      linewidth = 0.4,
      alpha = 0.1, legend = "none", outlier.shape = NA
    ) +
      geom_point(
        data = dat, aes(x = diagnosis_simple, y = !!sym(x), fill = diagnosis_simple, color = diagnosis_simple),
        size = 1.4,
        alpha = 0.5,
        position = position_jitter(w = 0.2)
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
      stat_pvalue_manual(
        pwc,
        label = "p.adj.str",
        tip.length = 0.01,
        label.size = 2.5,
        hide.ns = TRUE
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.1))
      ) +
      # # scale_x_discrete(labels = labels_xaxis, expand = c(0.5, 0.5)) +
      # ylim(min(dat[[sensor]]) - 0.1 * min(dat[[sensor]]), max(pwc$y.position) + 0.1 * max(pwc$y.position)) +
      labs(subtitle = paste0("Sensor: ", x)) +
      xlab("") +
      ylab("Intensity [a.u.]") +
      theme_classic() +
      ggtheme_no_legend
    bp
  }
)
names(bps) = sensors

plt_bps = plot_grid(
  plotlist = bps,
  ncol = 4,
  align = "hv"
)
plt_bps

