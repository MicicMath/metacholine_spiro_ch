source("./data_prep.R")

dat = test_set_final_ast_hc_cf_dat

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
      diff_ast_vs_hc_ci = paste0(
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

######
# test
######

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
temp_plt_s6bh = ggplot(dat, aes(x = diagnosis_simple, y = S6BH)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = jitter_pos) +
  geom_text_repel(
    aes(label = m_id),
    position = jitter_pos,
    max.overlap = Inf
  ) +
  labs(title = "S6BH") +
  theme_classic()

plot_grid(
  temp_plt_s3, temp_plt_s7
)
