source("./data_prep.R")

dat = final_ast_hc_cf_dat

###############################################################
# restrict to mid and high metacholine response only for asthma
###############################################################

table(dat$diagnosis_simple, dat$batch)

hc_tests = sapply(sensors, function(s) {
  wilcox.test(
    dat[[s]][dat$diagnosis_simple == "HC" & dat$batch == "batch_20260105"],
    dat[[s]][dat$diagnosis_simple == "HC" & dat$batch == "batch_20260808"]
  )$p.value
})

hc_dat = dat[dat$diagnosis_simple == "HC", ]

hc_p_adj = p.adjust(hc_tests, method = "fdr")
names(hc_p_adj) = sensors

hc_bps = lapply(
  sensors,
  function(x) {
    pwc = data.frame(
      group1 = "batch_20260105",
      group2 = "batch_20260808",
      y.position = max(hc_dat[[x]], na.rm = TRUE) * 1.1,
      p.adj = as.numeric(hc_p_adj[x]),
      p.adj.str = p_rnd(hc_p_adj[x])
    )
    bp = ggboxplot(hc_dat,
      x = "batch", y = x,
      width = 0.5,
      fill = "white",
      linewidth = 0.4,
      alpha = 0.1, legend = "none", outlier.shape = NA
    ) +
      geom_point(
        data = hc_dat, aes(x = batch, y = !!sym(x), fill = batch, color = batch),
        size = 2.4,
        alpha = 0.5,
        position = position_jitter(w = 0.2)
      ) +
      scale_color_manual(values = c("#228833", "#AA3377"), name = NULL) +
      scale_fill_manual(values = c("#228833", "#AA3377"), name = NULL) +
      stat_pvalue_manual(
        pwc,
        label = "p.adj.str",
        tip.length = 0.01,
        label.size = 3.5,
        hide.ns = TRUE
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      labs(subtitle = paste0("Sensor: ", x)) +
      xlab("") +
      ylab("Intensity [a.u.]") +
      theme_classic() +
      ggtheme_no_legend
    bp
  }
)
names(hc_bps) = sensors

plt_hc_bps = plot_grid(plotlist = hc_bps, ncol = 4, align = "hv")
plt_hc_bps

##############
# pca hc batch
##############

X_pca = prcomp(
  as.matrix(hc_dat[, sensors]),
  center = TRUE,
  scale. = TRUE
)

pca_importance = as.data.frame(summary(X_pca)$importance)

pca_dat = data.frame(
  X_pca$x,
  m_id = hc_dat$m_id,
  batch = hc_dat$batch
)

pc1_lab = paste0("PC1 ", round(pca_importance$PC1[2] * 100, 1), "%")
pc2_lab = paste0("PC2 ", round(pca_importance$PC2[2] * 100, 1), "%")
pc3_lab = paste0("PC3 ", round(pca_importance$PC3[2] * 100, 1), "%")
pc4_lab = paste0("PC4 ", round(pca_importance$PC4[2] * 100, 1), "%")

# diagnosis_colors = c(
#   HC  = "#4477AA",
#   AST = "#CC6677",
#   AST_AL = "#CC6677",
#   AST_NOAL = "#DDCC77"
# )

pca_by_batch_1 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC2, color = batch)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = batch), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by batch",
    x = pc1_lab,
    y = pc2_lab,
    color = "Batch"
  ) +
  theme_classic() +
  ggtheme_top_right

pca_by_batch_2 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC3, color = batch)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = batch), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by batch",
    x = pc1_lab,
    y = pc3_lab,
    color = "Batch"
  ) +
  theme_classic() +
  ggtheme_top_right

pca_by_batch_2b = ggplot(
  pca_dat,
  aes(x = PC2, y = PC3, color = batch)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = batch), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by batch",
    x = pc2_lab,
    y = pc3_lab,
    color = "Batch"
  ) +
  theme_classic() +
  ggtheme_top_right

pca_by_batch_3 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC4, color = batch)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = batch), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by batch",
    x = pc1_lab,
    y = pc4_lab,
    color = "Batch"
  ) +
  theme_classic() +
  ggtheme_top_right

plt_pca_hc = cowplot::plot_grid(
  pca_by_batch_1,
  pca_by_batch_2,
  pca_by_batch_2b,
  pca_by_batch_3,
  nrow = 2,
  align = "hv"
)
plt_pca_hc

#####################
# restrict to batch 2
#####################

# dat = dat[dat$batch == "batch_20260808", , drop = FALSE]

#########################
# remove S2 normalization
#########################

# X = as.matrix(dat[, sensors])
# X_clr = log(X) - rowMeans(log(X))
#
# dat[, sensors] = X_clr

####
# st
####

w_test = mclapply(
  sensors,
  function(x) {
    ast = dat[dat$diagnosis_simple == "AST", x]
    hc = dat[dat$diagnosis_simple == "HC", x]
    w = wilcox.test(ast, hc, conf.int = TRUE)
    cliff = cliff.delta(ast, hc, conf.level = 0.95)
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
      cliffs_delta_ci = paste0(
        n_rnd(cliff$estimate, 2),
        " (",
        n_rnd(cliff$conf.int[[1]], 2),
        ", ",
        n_rnd(cliff$conf.int[[2]], 2),
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
  "Est. difference (95% CI)",
  "Cliff's delta",
  "Raw p-value",
  "Adjusted p-value"
)

p_dat_format

# write.csv(p_dat_format, "./tables/univariate_results.csv", row.names = FALSE)

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
        size = 2.4,
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
        label.size = 3.5,
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

# ggsave("./plots/ast_hc_s3.png", bps$S3, width = 10, height = 10, units = "cm", dpi = 400, bg = "white")
# ggsave("./plots/ast_hc_s7.png", bps$S7, width = 10, height = 10, units = "cm", dpi = 400, bg = "white")

#############
# spider plot
#############

plot_sensor_pair = function(dat, sensor_x, sensor_y) {
  centers = aggregate(
    dat[c(sensor_x, sensor_y)],
    by = list(diagnosis_simple = dat$diagnosis_simple),
    FUN = mean
  )

  names(centers)[names(centers) == sensor_x] = paste0(sensor_x, "_center")
  names(centers)[names(centers) == sensor_y] = paste0(sensor_y, "_center")

  dat_plot = merge(
    dat,
    centers,
    by = "diagnosis_simple"
  )

  x_center = paste0(sensor_x, "_center")
  y_center = paste0(sensor_y, "_center")

  ggplot(
    dat_plot,
    aes(
      x = .data[[sensor_x]],
      y = .data[[sensor_y]],
      color = diagnosis_simple
    )
  ) +
    geom_segment(
      aes(
        x = .data[[x_center]],
        y = .data[[y_center]],
        xend = .data[[sensor_x]],
        yend = .data[[sensor_y]],
        color = diagnosis_simple
      ),
      alpha = 0.25,
      linewidth = 0.4
    ) +
    geom_point(size = 2.4, alpha = 0.65) +
    geom_point(
      aes(
        x = .data[[x_center]],
        y = .data[[y_center]],
        fill = diagnosis_simple
      ),
      shape = 21,
      size = 1,
      stroke = 0.5
    ) +
    scale_color_manual(
      values = c(AST = "#CC6677", HC = "#4477AA"),
      labels = c(AST = "Asthma", HC = "Healthy")
    ) +
    scale_fill_manual(
      values = c(AST = "#CC6677", HC = "#4477AA"),
      labels = c(AST = "Asthma", HC = "Healthy")
    ) +
    labs(
      x = paste0(sensor_x, " intensity [a.u.]"),
      y = paste0(sensor_y, " intensity [a.u.]")
    ) +
    theme_classic() +
    ggtheme_top_right
}

plt_s3_s7 = plot_sensor_pair(dat, "S7", "S3")
plt_s3_s7

# ggsave("./plots/spider_plot_s3_s7.png", plt_s3_s7, width = 10, height = 10, units = "cm", dpi = 400, bg = "white")

######
# test
######

# dat$S3_S7_logratio = log(dat$S3 / dat$S7)
#
# aggregate(
#   S3_S7_logratio ~ diagnosis_simple,
#   data = dat,
#   FUN = function(x) c(
#     mean = mean(x),
#     sd = sd(x),
#     IQR = IQR(x)
#   )
# )
#
# ggboxplot(
#   dat,
#   x = "diagnosis_simple",
#   y = "S3_S7_logratio",
#   add = "jitter"
# )
#
# intensity_sensors = c("S1", "S3", "S4", "S5", "S6", "S7")
#
# pairs = combn(intensity_sensors, 2, simplify = FALSE)
#
# ratio_variability = do.call(
#   rbind,
#   lapply(pairs, function(z) {
#
#     lr = log(dat[[z[1]]] / dat[[z[2]]])
#
#     data.frame(
#       pair = paste(z, collapse = "/"),
#       sd_HC = sd(lr[dat$diagnosis_simple == "HC"]),
#       sd_AST = sd(lr[dat$diagnosis_simple == "AST"])
#     )
#   })
# )
#
# ratio_variability

# plot_grid(
#   plt_hc_bps, plt_bps
# )

jitter_pos = position_jitter(
  width = 0.15,
  height = 0,
  seed = 123
)

temp_plt_s7 = ggplot(dat, aes(x = diagnosis_simple, y = S7)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = jitter_pos) +
  geom_text_repel(
    aes(label = m_id),
    position = jitter_pos,
    max.overlap = Inf
  ) +
  labs(title = "S7") +
  theme_classic()
temp_plt_s3 = ggplot(dat, aes(x = diagnosis_simple, y = S3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = jitter_pos) +
  geom_text_repel(
    aes(label = m_id),
    position = jitter_pos,
    max.overlap = Inf
  ) +
  labs(title = "S3") +
  theme_classic()

plot_grid(
  temp_plt_s7, temp_plt_s3
)
