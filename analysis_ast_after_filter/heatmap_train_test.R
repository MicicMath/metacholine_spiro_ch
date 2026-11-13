source("./data_prep.R")

train_dat = final_ast_hc_cf_dat
test_dat = test_set_final_ast_hc_cf_dat

################################
# within-dataset sensor z-scores
################################

# Each sensor is centred and scaled across the full training dataset and,
# separately, across the full test dataset, using AST and HC together.
# Colours describe relative levels within each dataset, not raw intensity.
# Diagnosis is used only for annotation, not scaling or clustering.

within_dataset_zscores = function(data) {
  X = as.matrix(data[, sensors, drop = FALSE])
  rownames(X) = as.character(data$m_id)

  if (any(!is.finite(X))) {
    stop("Heatmaps require finite sensor values.")
  }

  sensor_sd = apply(X, 2, sd)
  if (any(!is.finite(sensor_sd) | sensor_sd == 0)) {
    stop("Each sensor must have a finite, non-zero SD within the dataset.")
  }

  X_scaled = scale(X, center = TRUE, scale = TRUE)

  # Sensors in rows, measurements in columns.
  return(t(X_scaled))
}

z_train = within_dataset_zscores(train_dat)
z_test = within_dataset_zscores(test_dat)

##########################
# hierarchical clustering
##########################

# Euclidean distance and Ward linkage, as in simple_pca.R.
sensor_hclust = hclust(dist(z_train, method = "euclidean"), method = "ward.D2")
sensor_order = sensor_hclust$order

measurement_hclust_train = hclust(
  dist(t(z_train), method = "manhattan"),
  method = "complete"
)
measurement_hclust_test = hclust(
  dist(t(z_test), method = "manhattan"),
  method = "complete"
)

# Keep the training-derived sensor order in both plots for comparison.
# Test measurements are clustered independently; diagnosis does not
# determine their order or split them into predefined groups.

########################
# colours and annotations
########################

diagnosis_colors = c(AST = "#CC6677", HC = "#4477AA")

# Shared colour scale. Values outside +/-3 use the end colours;
# clustering still uses the full, unclipped z-scores.
heatmap_colors = circlize::colorRamp2(
  c(-3, 0, 3),
  c("#2166AC", "#FFFFFF", "#B2182B")
)

make_sensor_heatmap = function(data, z_mat, title, row_clustering, column_clustering) {
  stopifnot(identical(colnames(z_mat), as.character(data$m_id)))

  annotations = ComplexHeatmap::HeatmapAnnotation(
    Diagnosis = factor(data$diagnosis_simple, levels = c("AST", "HC")),
    col = list(Diagnosis = diagnosis_colors),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontsize = 9),
    simple_anno_size = grid::unit(3, "mm"),
    annotation_legend_param = list(
      Diagnosis = list(at = c("AST", "HC"), labels = c("Asthma", "Healthy controls"))
    )
  )

  ComplexHeatmap::Heatmap(
    z_mat,
    name = "Sensor z-score",
    col = heatmap_colors,
    top_annotation = annotations,
    # cluster_rows = row_clustering,
    row_order = sensor_order,
    row_dend_reorder = FALSE,
    cluster_columns = column_clustering,
    column_dend_reorder = FALSE,
    row_dend_width = grid::unit(12, "mm"),
    column_dend_height = grid::unit(15, "mm"),
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 10),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Sensor z-score\n(within dataset)",
      at = c(-3, 0, 3),
      labels = c("<= -3", "0", ">= 3")
    ),
    use_raster = FALSE
  )
}

heatmap_train = make_sensor_heatmap(
  data = train_dat,
  z_mat = z_train,
  title = paste0("Training (n = ", nrow(train_dat), ")"),
  row_clustering = sensor_hclust,
  column_clustering = measurement_hclust_train
)

# The test heatmap reuses the training row order, without displaying a
# sensor dendrogram that could be mistaken for a test-derived clustering.
heatmap_test = make_sensor_heatmap(
  data = test_dat,
  z_mat = z_test,
  title = paste0("Test (n = ", nrow(test_dat), "; sensor order from training)"),
  row_clustering = FALSE,
  column_clustering = measurement_hclust_test
)

####################
# convert to ggplot
####################

heatmap_as_ggplot = function(heatmap) {
  # Capture at a landscape size so annotations and legends have room.
  heatmap_grob = grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heatmap,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE
    ),
    width = 10,
    height = 5
  )

  ggplotify::as.ggplot(heatmap_grob)
}

plt_heatmap_train = heatmap_as_ggplot(heatmap_train)
plt_heatmap_test = heatmap_as_ggplot(heatmap_test)

plt_heatmaps = cowplot::plot_grid(
  plt_heatmap_train,
  plt_heatmap_test,
  ncol = 1
)
print(plt_heatmaps)

############
# save plots
############

# ggsave(
#   filename = "./plots/heatmap_train_test.png",
#   plot = plt_heatmaps,
#   width = 28,
#   height = 28,
#   units = "cm",
#   dpi = 300,
#   bg = "white"
# )
