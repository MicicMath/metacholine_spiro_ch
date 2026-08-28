source("./data_prep.R")

dat = final_ast_hc_cf_dat

#################
# subset to group
#################

# dat = dat[dat$diagnosis_simple == "AST", , drop = FALSE]

#####
# cor
#####

M = cor(dat[, sensors], use = "pairwise.complete.obs")

cor_plot = function(data, title) {
  M = cor(data, use = "pairwise.complete.obs")
  hc = hclust(dist(1 - M), method = "average")
  ord = hc$order

  M_long = reshape2::melt(M)

  M_long = M_long[
    as.numeric(M_long$Var1) < as.numeric(M_long$Var2),
  ]

  plt_cor = ggplot(M_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white", linewidth = 1) +
    
    geom_text(
      aes(label = sprintf("%.2f", value)),
      size = 3,
      color = "black"
    ) +
    
    scale_fill_gradient2(
      low = "#BB4444",
      mid = "#FFFFFF",
      high = "#4477AA",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    
    # coord_fixed() +
    scale_x_discrete(position = "top") +

    labs(title = title) +
    
    theme_minimal(base_size = 12) +
    theme(
      # legend.position = "bottom",
      axis.title = element_blank(),
      panel.grid = element_blank(),
      # axis.text.x = element_text(angle = 45, hjust = .5, vjust = 1),
      axis.text.x = element_text(),
      axis.text.y = element_text()
    )

  plt_cor
}

plt_cor_1 = cor_plot(dat[dat$diagnosis_simple == "AST", sensors], "Asthma group")
plt_cor_2 = cor_plot(dat[dat$diagnosis_simple == "HC", sensors], "Healthy controls")

plot_grid(
  plt_cor_1,
  plt_cor_2
)

#############################################
# permutation test between asthma and healthy
#############################################

M_ast = cor(
  dat[dat$diagnosis_simple == "AST", sensors],
  use = "pairwise.complete.obs"
)

M_hc = cor(
  dat[dat$diagnosis_simple == "HC", sensors],
  use = "pairwise.complete.obs"
)

# observed global difference between correlation matrices
idx = upper.tri(M_ast)

T_obs = sum((M_ast[idx] - M_hc[idx])^2)

# permutation test
set.seed(101)

B = 10000

diagnosis = dat$diagnosis_simple

T_perm = numeric(B)

for (b in 1:B) {

  diagnosis_perm = sample(diagnosis)

  M_ast_perm = cor(
    dat[diagnosis_perm == "AST", sensors],
    use = "pairwise.complete.obs"
  )

  M_hc_perm = cor(
    dat[diagnosis_perm == "HC", sensors],
    use = "pairwise.complete.obs"
  )

  T_perm[b] = sum(
    (M_ast_perm[idx] - M_hc_perm[idx])^2
  )
}

# permutation p-value
p_value = (sum(T_perm >= T_obs) + 1) / (B + 1)

T_obs
p_value

#####
# pca
#####

pca_plot = function(data, comp1 = "PC1", comp2 = "PC2", title) {
  X_pca = prcomp(
    as.matrix(data[, sensors]),
    center = TRUE,
    scale. = TRUE
  )

  pca_importance = as.data.frame(summary(X_pca)$importance)

  pca_dat = data.frame(
    X_pca$x,
    m_id = data$m_id
  )

  comp1_lab = paste0(comp1, " ", round(pca_importance[[comp1]][2] * 100, 1), "%")
  comp2_lab = paste0(comp2, " ", round(pca_importance[[comp2]][2] * 100, 1), "%")

  plt = ggplot(
    pca_dat,
    aes(x = !!sym(comp1), y = !!sym(comp2))
  ) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
    # stat_ellipse(aes(group = batch), show.legend = FALSE, linewidth = 0.5) +
    geom_point(size = 2.4, alpha = 0.75) +
    # scale_color_manual(values = diagnosis_colors) +
    labs(
      title = paste0("PCA scores ", title),
      x = comp1_lab,
      y = comp2_lab
    ) +
    theme_classic() +
    ggtheme_top_right

    return(plt)
}

pc_ast_1 = pca_plot(dat[dat$diagnosis_simple == "AST", ], "PC1", "PC2", "asthma group")
pc_hc_1 = pca_plot(dat[dat$diagnosis_simple == "HC", ], "PC1", "PC2", "healty control")

plot_grid(
  pc_hc_1,
  pc_ast_1
)
