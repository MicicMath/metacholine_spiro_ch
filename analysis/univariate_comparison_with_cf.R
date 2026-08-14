source("./data_prep.R")

dat = final_ast_hc_cf_dat[, cols]

dat$diagnosis = ifelse(
  dat$diagnosis == "PUL_ASTHMA_ALLERGIC",
  "AST_AL",
  ifelse(dat$diagnosis == "PUL_ASTHMA_NONALLERGIC",
    "AST_NOAL",
    ifelse(dat$diagnosis == "PUL_CF",
      "CF", "HC"
    )
  )
)

#############################
# exclude asthma non-allergic
#############################

dat = dat[dat$diagnosis != "AST_NOAL", , drop = FALSE]

#######################
# remap diagnosis label
#######################

dat$diagnosis = ifelse(
  grepl("AST", dat$diagnosis),
  "AST", dat$diagnosis
)

#######################
# add metacholine label
#######################

dat$metacho = ifelse(
  (dat$diagnosis == "AST") & (dat$metacholine_test_result == 0),
  "R0", ifelse(
    (dat$diagnosis == "AST") & (dat$metacholine_test_result == 1),
    "R1", ifelse(
      (dat$diagnosis == "AST") & (dat$metacholine_test_result == 2),
      "R2", ifelse(
        (dat$diagnosis == "AST") & (dat$metacholine_test_result == 3),
        "R3", dat$diagnosis
      )
    )
  )
)

###################################
# remove asthma with no metacholine
###################################

mask_no_metacho = (dat$diagnosis == "AST") & is.na(dat$metacholine_test_result)

dat = dat[!mask_no_metacho, , drop = FALSE]

#################################
# restrict to confirmed diagnosis
#################################

conf_diag = (dat$diagnosis_status == 1) | is.na(dat$diagnosis_status)
dat = dat[conf_diag, , drop = FALSE]

#######################################
# restrict to metacholine level 2 and 3
#######################################

metacho_23 = (dat$metacho %in% c("R2", "R3")) | dat$diagnosis == "HC" | dat$diagnosis == "CF"
dat = dat[metacho_23, , drop = FALSE]

####
# st
####

post_hoc = mclapply(
  sensors,
  function(x) {
    form = as.formula(paste0(x, " ~ diagnosis"))
    krusk = kruskal_test(data = dat, form)

    dt = dunn_test(form, data = dat, p.adjust.method = "holm", detailed = TRUE)

    pwc = pairwise_wilcox_test(dat, form, p.adjust.method = "bonferroni", detailed = TRUE)
    factor_levels = levels(dat$diagnosis)
    pwc = pwc[order(match(pwc$group1, factor_levels), match(pwc$group2, factor_levels)), ]

    max_y = max(dat[[x]], na.rm = TRUE)
    pwc$y.position = max_y + seq(0.1, 0.5, length.out = nrow(pwc)) * max_y

    expected_W = (pwc$n1 * pwc$n2) / 2 # Expected W under null
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

color_alt = c(
  "#CC6677",
  "#AA4499",
  "#4477AA"
)

bp_plots = function(res, dat, sensor, p_adj) {
  pd = ggplot2::position_jitter(w = 0.1, seed = 1)
  name = sensor
  deg = round(res$krusk$statistic, 2)
  num = res$krusk$n

  p_val = p_rnd(res$krusk$p)
  p_adj = p_rnd(p_adj[names(p_adj) == sensor])

  pw_dunn = res$dunn
  pw_dunn$p.adj = sapply(pw_dunn$p.adj, function(x) p_rnd(x))

  pw_dunn$y.position = max(dat[[sensor]], na.rm = TRUE) +
    seq(0.1, 0.5, length.out = nrow(pw_dunn)) * max(dat[[sensor]], na.rm = TRUE)

  bp1 = ggboxplot(dat,
    x = "diagnosis", y = sensor,
    width = 0.5,
    fill = "lightgrey",
    color = "grey",
    alpha = 0.1, legend = "none", outlier.shape = NA
  ) +
    geom_point(
      data = dat, aes(x = diagnosis, y = !!sym(sensor), fill = diagnosis, color = diagnosis),
      size = 1, alpha = 0.5, position = position_jitter(w = 0.2)
    ) +
    scale_color_manual(values = color_alt) +
    scale_fill_manual(values = color_alt) +
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
      # title = name,
      subtitle = substitute(
        paste(
          # "Kruskal-Wallis, ",
          italic("p")["adj"], ": ", p_adj
        ),
        list(
          p_val = p_val,
          p_adj = p_adj
        )
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


