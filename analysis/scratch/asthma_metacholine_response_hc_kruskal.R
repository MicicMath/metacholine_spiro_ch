source("./data_prep.R")

dat = final_ast_hc_cf_dat
dat = dat[dat$diagnosis_simple %in% c("AST", "HC"), , drop = FALSE]

#########################################################
# compare HC, negative AST, and metacholine-response AST
#########################################################

response_mask = dat$diagnosis_simple == "AST" & dat$metacholine_response %in% c("Mid", "High")
neg_mask = dat$diagnosis_simple == "AST" & dat$metacholine_response == "Neg"
hc_mask = dat$diagnosis_simple == "HC"

dat = dat[response_mask | neg_mask | hc_mask, , drop = FALSE]

dat$metacholine_hc_group = ifelse(
  dat$diagnosis_simple == "HC",
  "HC",
  ifelse(dat$metacholine_response == "Neg", "AST_Neg", "AST_Response")
)
dat$metacholine_hc_group = factor(
  dat$metacholine_hc_group,
  levels = c("HC", "AST_Neg", "AST_Response")
)

table(dat$metacholine_hc_group)

####
# st
####

post_hoc = mclapply(
  sensors,
  function(x) {
    form = as.formula(paste0(x, " ~ metacholine_hc_group"))
    krusk = kruskal_test(data = dat, form)

    dt = dunn_test(form, data = dat, p.adjust.method = "holm", detailed = TRUE)

    pwc = pairwise_wilcox_test(dat, form, p.adjust.method = "bonferroni", detailed = TRUE)
    factor_levels = levels(dat$metacholine_hc_group)
    pwc = pwc[order(match(pwc$group1, factor_levels), match(pwc$group2, factor_levels)), ]

    max_y = max(dat[[x]], na.rm = TRUE)
    pwc$y.position = max_y + seq(0.1, 0.5, length.out = nrow(pwc)) * max_y

    expected_W = (pwc$n1 * pwc$n2) / 2
    sd_W = sqrt((pwc$n1 * pwc$n2 * (pwc$n1 + pwc$n2 + 1)) / 12)
    Z = (pwc$statistic - expected_W) / sd_W
    N = pwc$n1 + pwc$n2
    effs = abs(Z) / sqrt(N)

    pwc_dat = data.frame(
      feat_id = x,
      comparison = paste0(pwc$group1, " - ", pwc$group2),
      estimate = paste0(round(pwc$estimate, 2), " (", round(pwc$conf.low, 2), ", ", round(pwc$conf.high, 2), ")"),
      p_post_hoc = sapply(pwc$p, function(x) p_rnd(x)),
      p_adj_post_hoc = sapply(pwc$p.adj, function(x) p_rnd(x)),
      sig = pwc$p.adj.signif,
      eff = sapply(effs, function(x) round(x, 3)),
      eff_intr = sapply(effs, function(x) eff_interpretation(x))
    )

    list(
      krusk = data.frame(krusk),
      dunn = data.frame(dt),
      wlcx = pwc,
      pwc_dat = pwc_dat
    )
  }, mc.cores = parallel::detectCores() - 1
)
names(post_hoc) = sensors

p_krusk = sapply(post_hoc, function(x) x$krusk$p)
p_adj_krusk = p.adjust(p_krusk, method = "fdr")

krusk_dat_format = data.frame(
  Sensor = sensors,
  p_val = sapply(p_krusk, function(x) p_rnd(x)),
  p_adj = sapply(p_adj_krusk, function(x) p_rnd(x))
)
krusk_dat_format = krusk_dat_format[order(p_adj_krusk), , drop = FALSE]
colnames(krusk_dat_format) = c("Sensor", "Raw p-value", "Adjusted p-value")
krusk_dat_format

dunn_post_hoc = lapply(
  post_hoc,
  function(x) {
    x$dunn
  }
)
dunn_post_hoc = do.call(rbind, dunn_post_hoc)
rownames(dunn_post_hoc) = NULL

#####
# bps
#####

bp_plots = function(res, dat, sensor, p_adj) {
  p_adj = p_rnd(p_adj[names(p_adj) == sensor])

  pw_dunn = res$dunn
  pw_dunn$p.adj = sapply(pw_dunn$p.adj, function(x) p_rnd(x))

  pw_dunn$y.position = max(dat[[sensor]], na.rm = TRUE) +
    seq(0.1, 0.5, length.out = nrow(pw_dunn)) * max(dat[[sensor]], na.rm = TRUE)

  bp1 = ggboxplot(dat,
    x = "metacholine_hc_group", y = sensor,
    width = 0.5,
    fill = "lightgrey",
    color = "grey",
    alpha = 0.1, legend = "none", outlier.shape = NA
  ) +
    geom_point(
      data = dat, aes(x = metacholine_hc_group, y = !!sym(sensor), fill = metacholine_hc_group, color = metacholine_hc_group),
      size = 1, alpha = 0.5, position = position_jitter(w = 0.2)
    ) +
    scale_color_manual(values = c(HC = "#4477AA", AST_Neg = "#999999", AST_Response = "#CC6677")) +
    scale_fill_manual(values = c(HC = "#4477AA", AST_Neg = "#999999", AST_Response = "#CC6677")) +
    stat_pvalue_manual(
      pw_dunn,
      label = "p.adj",
      tip.length = 0,
      label.size = 2.5,
      hide.ns = TRUE
    ) +
    ylim(min(dat[[sensor]]) - 0.1 * min(dat[[sensor]]), max(pw_dunn$y.position) + 0.1 * max(pw_dunn$y.position)) +
    xlab("") +
    labs(
      subtitle = substitute(
        paste(italic("p")["adj"], ": ", p_adj),
        list(p_adj = p_adj)
      )
    ) +
    theme_classic() +
    ggtheme_no_legend_45
  return(bp1)
}

plts = lapply(
  sensors,
  function(x) {
    bp_plots(post_hoc[[x]], dat, x, p_adj_krusk)
  }
)
names(plts) = sensors
plts = plts[order(p_adj_krusk)]

plt_all = plot_grid(
  plotlist = plts,
  ncol = 7,
  align = "hv"
)
plt_all

# ggsave("./plots/ast_metacholine_response_hc_kruskal.png", plt_all, width = 24, height = 12, units = "cm", dpi = 400, bg = "white")
