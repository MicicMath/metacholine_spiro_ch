source("./data_prep.R")

train_dat = final_ast_hc_cf_dat
test_dat = test_set_final_ast_hc_cf_dat

X_train = train_dat[, sensors]
X_test = test_dat[, sensors]

meta_data_train = train_dat[, setdiff(colnames(train_dat), sensors)]
meta_data_test = test_dat[, setdiff(colnames(test_dat), sensors)]

diagnosis_colors = c(AST = "#CC6677", HC = "#4477AA")

#####
# pca
#####

run_pca_dat = function(X_mat, meta_data) {
  var_cols = apply(X_mat, 2, var, na.rm = TRUE)
  valid = !is.na(var_cols) & var_cols > 0
  pca = prcomp(X_mat[, valid], center = TRUE, scale. = TRUE)
  var_exp = (pca$sdev^2) / sum(pca$sdev^2) * 100

  dat_pca = data.frame(
    pca$x,
    label = meta_data$diagnosis_simple,
    patient_id = meta_data$patient_id
  )
  list(dat = dat_pca, var = var_exp)
}

pca_train = run_pca_dat(X_train, meta_data_train)

plt_pca_train = ggplot(pca_train$dat, aes(x = PC1, y = PC2, color = label)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  # geom_text_repel(aes(label = patient_id)) +
  scale_color_manual(values = diagnosis_colors, name = NULL) +
  scale_fill_manual(values = diagnosis_colors, name = NULL) +
  labs(
    title = "PCA train data",
    x = sprintf("PC1 (%.1f%%)", pca_train$var[1]),
    y = sprintf("PC2 (%.1f%%)", pca_train$var[2])
  ) +
  theme_bw()
plt_pca_train

pca_test = run_pca_dat(X_test, meta_data_test)

plt_pca_test = ggplot(pca_test$dat, aes(x = PC1, y = PC2, color = label)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  # geom_text_repel(aes(label = patient_id)) +
  scale_color_manual(values = diagnosis_colors, name = NULL) +
  scale_fill_manual(values = diagnosis_colors, name = NULL) +
  labs(
    title = "PCA test data",
    x = sprintf("PC1 (%.1f%%)", pca_test$var[1]),
    y = sprintf("PC2 (%.1f%%)", pca_test$var[2])
  ) +
  theme_bw()
plt_pca_test

plot_grid(
  plt_pca_train, plt_pca_test
)
