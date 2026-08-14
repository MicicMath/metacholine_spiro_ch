source("./data_prep.R")

dat = final_ast_hc_dat[, cols]

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

########
# counts
########

table(dat$diagnosis, dat$visit_exhal)
table(dat$diagnosis, dat$metacholine_test_result)

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
  metacho = dat$metacho,
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
  aes(x = PC1, y = PC2, color = metacho)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = metacho), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by diagnosis",
    x = pc1_lab,
    y = pc2_lab
    # color = "Diagnosis"
  ) +
  theme_classic()

pca_by_diagnosis_2 = ggplot(
  pca_dat,
  aes(x = PC1, y = PC3, color = metacho)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70", linetype = "dashed") +
  stat_ellipse(aes(group = metacho), show.legend = FALSE, linewidth = 0.5) +
  geom_point(size = 2.4, alpha = 0.75) +
  # scale_color_manual(values = diagnosis_colors) +
  labs(
    title = "PCA scores by diagnosis",
    x = pc1_lab,
    y = pc3_lab
    # color = "Diagnosis"
  ) +
  theme_classic()

plt_pca = cowplot::plot_grid(
  pca_by_diagnosis_1,
  pca_by_diagnosis_2,
  nrow = 1,
  align = "hv"
)
plt_pca
