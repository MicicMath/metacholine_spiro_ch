source("./data_prep.R")

dat = final_ast_hc_dat[, cols]

dat$diagnosis = ifelse(
  dat$diagnosis == "PUL_ASTHMA_ALLERGIC",
  "AST_AL",
  ifelse(dat$diagnosis == "PUL_ASTHMA_NONALLERGIC",
    "AST_NOAL", "HC"
  )
)

dat$diagnosis_simple = ifelse(dat$diagnosis == "HC", "HC", "AST")

#################################
# restrict to confirmed diagnosis
#################################

conf_diag = (dat$diagnosis_status == 1) | is.na(dat$diagnosis_status)

#####
# pca
#####

X_pca = prcomp(
  as.matrix(dat[, sensors]),
  center = TRUE,
  scale. = TRUE
)

pca_importance = as.data.frame(summary(X_pca)$importance)

pca_dat = data.frame(
  X_pca$x,
  m_id = dat$m_id,
  diagnosis = dat$diagnosis
)

pc1_lab = paste0("PC1 ", round(pca_importance$PC1[2] * 100, 1), "%")
pc2_lab = paste0("PC2 ", round(pca_importance$PC2[2] * 100, 1), "%")
pc3_lab = paste0("PC3 ", round(pca_importance$PC3[2] * 100, 1), "%")

diagnosis_colors = c(
  HC  = "#4477AA",
  AST_AL = "#CC6677",
  AST_NOAL = "#DDCC77"
)

pca_by_diagnosis_1 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC2, color = diagnosis)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = diagnosis), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by diagnosis",
    x = pc1_lab,
    y = pc2_lab,
    color = "Diagnosis"
  ) +
  theme_classic()

pca_by_diagnosis_2 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC3, color = diagnosis)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = diagnosis), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by diagnosis",
    x = pc1_lab,
    y = pc3_lab,
    color = "Diagnosis"
  ) +
  theme_classic()

plt_pca = cowplot::plot_grid(
  pca_by_diagnosis_1,
  pca_by_diagnosis_2,
  nrow = 1,
  align = "hv"
)
plt_pca
