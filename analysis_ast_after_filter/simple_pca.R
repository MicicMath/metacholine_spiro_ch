source("./data_prep.R")

dat = final_ast_hc_cf_dat

#######
# setup
#######

run_pca_dat = function(data) {
  X_mat = data[, sensors]
  var_cols = apply(X_mat, 2, var, na.rm = TRUE)
  valid = !is.na(var_cols) & var_cols > 0
  pca = prcomp(X_mat[, valid], center = TRUE, scale. = TRUE)
  var_exp = (pca$sdev^2) / sum(pca$sdev^2) * 100

  dat_pca = data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    label = data$diagnosis_simple,
    sample_id = data$m_id
  )
  list(d = dat_pca, var = var_exp)
}

#####
# pca
#####

pca_dat = run_pca_dat(data = dat)

p1 = ggplot(pca_dat$d, aes(x = PC1, y = PC2, color = label)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_ellipse(level = 0.95, linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(values = c("#CC6677", "#4477AA"), labels = c("Asthma", "Healthy controls")) +
  labs(
    title = "",
    x = sprintf("PC1 (%.1f%%)", pca_dat$var[1]),
    y = sprintf("PC2 (%.1f%%)", pca_dat$var[2])
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(1, 1),
    legend.justification = c("right", "top"),
    legend.key.spacing.y = unit(-0.2, "cm")
  )

p2 = ggplot(pca_dat$d, aes(x = PC1, y = PC3, color = label)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  stat_ellipse(level = 0.95, linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(values = c("#CC6677", "#4477AA"), labels = c("Asthma", "Healthy controls")) +
  labs(
    title = "",
    x = sprintf("PC1 (%.1f%%)", pca_dat$var[1]),
    y = sprintf("PC3 (%.1f%%)", pca_dat$var[3])
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(1, 1),
    legend.justification = c("right", "top"),
    legend.key.spacing.y = unit(-0.2, "cm")
  )

plot_grid(
  p1, p2
)

#####
# hierarchical clustering
#####

run_hclust_dat = function(data) {
  X_mat = data[, sensors]

  # Remove zero-variance sensors
  var_cols = apply(X_mat, 2, var, na.rm = TRUE)
  valid = !is.na(var_cols) & var_cols > 0
  X_mat = X_mat[, valid]

  # Standardize sensors
  X_scaled = scale(X_mat)

  # Hierarchical clustering
  d = dist(X_scaled, method = "euclidean")
  hc = hclust(d, method = "ward.D2")

  list(
    hc = hc,
    X = X_scaled
  )
}

hclust_dat = run_hclust_dat(data = dat)

p3 = factoextra::fviz_dend(
  hclust_dat$hc,
  k = 2,
  k_colors = c("#CC6677", "#4477AA"),
  rect = TRUE,
  rect_fill = FALSE,
  show_labels = FALSE,
  color_labels_by_k = FALSE
) +
  labs(
    title = "",
    x = "Samples",
    y = "Height"
  ) +
  theme_classic()
p3
