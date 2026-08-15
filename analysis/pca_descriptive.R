source("./data_prep.R")

dat = final_ast_hc_cf_dat
dat = dat[dat$diagnosis != "CF", , drop = FALSE]

###############################################################
# restrict to mid and high metacholine response only for asthma
###############################################################

mid_high_metach_response_mask = dat$metacholine_response %in% c("Mid", "High")
other_mask = dat$diagnosis_simple != "AST"

dat = dat[mid_high_metach_response_mask | other_mask, , drop = FALSE]

########
# counts
########

table(dat$diagnosis_simple, dat$visit_exhal)
table(dat$diagnosis_simple, dat$metacholine_response)

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
  metacholine_response = dat$metacholine_response,
  diagnosis = dat$diagnosis_simple
)

pc1_lab = paste0("PC1 ", round(pca_importance$PC1[2] * 100, 1), "%")
pc2_lab = paste0("PC2 ", round(pca_importance$PC2[2] * 100, 1), "%")
pc3_lab = paste0("PC3 ", round(pca_importance$PC3[2] * 100, 1), "%")

# diagnosis_colors = c(
#   HC  = "#4477AA",
#   AST = "#CC6677",
#   AST_AL = "#CC6677",
#   AST_NOAL = "#DDCC77"
# )

pca_by_diagnosis_1 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC2, color = metacholine_response)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = metacholine_response), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by diagnosis",
    x = pc1_lab,
    y = pc2_lab,
    color = "Metacholine\nresponse"
  ) +
  theme_classic()

pca_by_diagnosis_2 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC3, color = metacholine_response)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = metacholine_response), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by diagnosis",
    x = pc1_lab,
    y = pc3_lab,
    color = "Metacholine\nresponse"
  ) +
  theme_classic()

plt_pca = cowplot::plot_grid(
  pca_by_diagnosis_1,
  pca_by_diagnosis_2,
  nrow = 1,
  align = "hv"
)
plt_pca
