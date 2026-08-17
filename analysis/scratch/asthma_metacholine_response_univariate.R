source("./data_prep.R")

dat = final_ast_hc_cf_dat
dat = dat[dat$diagnosis_simple == "AST", , drop = FALSE]

#######################################################
# compare negative and positive metacholine response AST
#######################################################

response_mask = dat$metacholine_response %in% c("Mid", "High")
neg_mask = dat$metacholine_response == "Neg"

dat = dat[response_mask | neg_mask, , drop = FALSE]

dat$metacholine_group = ifelse(
  dat$metacholine_response == "Neg",
  "Neg",
  "Response"
)
dat$metacholine_group = factor(dat$metacholine_group, levels = c("Neg", "Response"))

table(dat$metacholine_group)
table(dat$metacholine_response)

####
# st
####

w_test = mclapply(
  sensors,
  function(x) {
    neg = dat[dat$metacholine_group == "Neg", x]
    response = dat[dat$metacholine_group == "Response", x]
    w = wilcox.test(response, neg, conf.int = TRUE)
    d = data.frame(
      sensor = x,
      diff_response_vs_neg_ci = paste0(
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
  },
  mc.cores = parallel::detectCores() - 2
)
w_test = do.call(rbind, w_test)

p_adj = p.adjust(w_test$p_val, method = "fdr")
names(p_adj) = sensors
w_test$p_adj = p_adj

p_dat_format = w_test
p_dat_format$p_val = p_rnd(p_dat_format$p_val)
p_dat_format$p_adj = p_rnd(p_dat_format$p_adj)

colnames(p_dat_format) = c(
  "Sensor",
  "Est. difference Response - Neg (95% CI)",
  "Raw p-value",
  "Adjusted p-value"
)

p_dat_format

####
# bp
####

bps = lapply(
  sensors,
  function(x) {
    pwc = data.frame(
      group1 = "Neg",
      group2 = "Response",
      y.position = max(dat[[x]], na.rm = TRUE) * 1.1,
      p.adj = as.numeric(p_adj[x]),
      p.adj.str = p_rnd(p_adj[x])
    )

    bp = ggboxplot(dat,
      x = "metacholine_group", y = x,
      width = 0.5,
      fill = "white",
      linewidth = 0.4,
      alpha = 0.1, legend = "none", outlier.shape = NA
    ) +
      geom_point(
        data = dat, aes(x = metacholine_group, y = !!sym(x), fill = metacholine_group, color = metacholine_group),
        size = 2.4,
        alpha = 0.5,
        position = position_jitter(w = 0.2)
      ) +
      scale_color_manual(
        values = c("#4477AA", "#CC6677"),
        name = NULL
      ) +
      scale_fill_manual(
        values = c("#4477AA", "#CC6677"),
        name = NULL
      ) +
      stat_pvalue_manual(
        pwc,
        label = "p.adj.str",
        tip.length = 0.01,
        label.size = 3.5,
        hide.ns = TRUE
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.1))
      ) +
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

# ggsave("./plots/ast_metacholine_response_univariate.png", plt_bps, width = 24, height = 18, units = "cm", dpi = 400, bg = "white")
